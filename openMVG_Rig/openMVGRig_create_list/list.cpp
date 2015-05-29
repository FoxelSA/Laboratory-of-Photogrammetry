/*
 * listomvg
 *
 * Copyright (c) 2014-2015 FOXEL SA - http://foxel.ch
 * Please read <http://foxel.ch/license> for more information.
 *
 *
 * Author(s):
 *
 *      Stéphane Flotron <s.flotron@foxel.ch>
 *
 * This file is part of the FOXEL project <http://foxel.ch>.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Additional Terms:
 *
 *      You are required to preserve legal notices and author attributions in
 *      that material or in the Appropriate Legal Notices displayed by works
 *      containing it.
 *
 *      You are required to attribute the work as explained in the "Usage and
 *      Attribution" section of <http://foxel.ch/license>.
 */

 /*! \file list.cpp
 * \author Stephane Flotron <s.flotron@foxel.ch>
 */
 /*! \mainpage listomvg
 * \section listomvg
 *
 * For Elphel's camera, generate file lists.txt needed by openMVG, for rigid rig and standard openMVG.
 *
 * \section Documentation
 *
 * Documentation can be consulted on the [wiki](https://github.com/FoxelSA/listomvg/wiki).
 *
 * \section Copyright
 *
 * Copyright (c) 2014-2015 FOXEL SA - [http://foxel.ch](http://foxel.ch)<br />
 * This program is part of the FOXEL project <[http://foxel.ch](http://foxel.ch)>.
 *
 * Please read the [COPYRIGHT.md](COPYRIGHT.md) file for more information.
 *
 * \section License
 *
 * This program is licensed under the terms of the
 * [GNU Affero General Public License v3](http://www.gnu.org/licenses/agpl.html)
 * (GNU AGPL), with two additional terms. The content is licensed under the terms
 * of the [Creative Commons Attribution-ShareAlike 4.0 International](http://creativecommons.org/licenses/by-sa/4.0/)
 * (CC BY-SA) license.
 *
 * You must read <[http://foxel.ch/license](http://foxel.ch/license)> for more
 *information about our Licensing terms and our Usage and Attribution guidelines.
 *
 */

#include <tools.hpp>
#include <types.hpp>
#include <list_utils.hpp>
#include <stdlib.h>
#include <ctype.h>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <stlplus3/filesystemSimplified/file_system.hpp>
#include <cmdLine/cmdLine.h>

using namespace std;

/*! \brief Main software function
*
* This function takes a sensor as input and load all calibration
* needed for gnomonic projection.
*
* \param sImageDir      Path/directory containing images
* \param smacAddress    Mac address of the elphel camera that take the photo
* \param sOutputDir     Path/directory where you want the file lists.txt lies
* \param sMountPoint    mount point of the camera folder
* \param sChannelFile   Optionnal Argument. Path/name of the channel file containing the subcamera you want to keep
*                       By default we keep all subcameras.
* \param bRigidRig      Optionnal Argument. if 0, do not use rigid rig structure, if 1
*                       use rigid rig structure. Default is 1.
* \param focalPixPermm  Optionnal Argument. If passed to the programm, we will use the
*                       focal given in input for all images. Warning, the focal should be in
*                       pixels per mm.
* \param bUseCalibPrincipalPoint Optionnal Argument. If 0, use the center of as principal point, if 1
*                       use the principal point given by calibration process.
*
* \return 0 if all was well, 1 in other cases.
*/


int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDir,
    sChannelFile = "",
    smacAddress = "",
    sOutputDir = "",
    sMountPoint = "",
    sTimestampLow= "",
    sTimestampUp="",
    sGpsFileName="";

  li_Real_t     focalPixPermm = -1.0;
  bool          bRigidRig     = true;
  bool          bUseCalibPrincipalPoint = false;

  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('m', smacAddress, "macAddress") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('d', sMountPoint, "mountPoint") );
  cmd.add( make_option('c', sChannelFile, "channelFile") );
  cmd.add( make_option('r', bRigidRig, "rigidRig") );
  cmd.add( make_option('p', bUseCalibPrincipalPoint, "useCalibPrincipalPoint") );
  cmd.add( make_option('f', focalPixPermm, "focal") );
  cmd.add( make_option('a', sTimestampLow, "lowerBound") );
  cmd.add( make_option('b', sTimestampUp, "upperBound") );
  cmd.add( make_option('g', sGpsFileName, "gps"));

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-m|--macAddress\n"
      << "[-o|--outputDirectory]\n"
      << "[-d|--mountPoint]\n"
      << "[-c|--channelFile]\n"
      << "[-r|--rigidRig \n"
      << "   -r 0 : no rigid rig \n"
      << "   -r 1 : with rigid rig structure\n"
      << "[-p|--useCalibPrincipalPoint\n"
      << "   -p 0 : do not use calibration principal point \n"
      << "   -p 1 : use calibration principal point \n"
      << "[-a|--lowerBound \n"
      << "[-b|--upperBound \n"
      << "[-f|--focal] (pixels)\n"
      << "[-g|--gps] GPU/IMU json file\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  // verify that input argumet are valid
  const bool bValidInput = isInputValid
                              (
                                    argv[0],
                                    sImageDir,
                                    sOutputDir,
                                    smacAddress,
                                    sMountPoint,
                                    sChannelFile,
                                    bRigidRig,
                                    bUseCalibPrincipalPoint,
                                    focalPixPermm,
                                    sTimestampLow,
                                    sTimestampUp
                              );
   if( !bValidInput )
   {
       std::cerr << " Input data are not valid. Please check your input\n" << std::endl;
       return EXIT_FAILURE;
   }
   else  // if input data is valid, go ahead
   {

       // now extract calibration information related to each module
       std::vector < sensorData > vec_sensorData;
       const bool bLoadCalibration = loadCalibrationData( vec_sensorData, sMountPoint, smacAddress);

       if( !bLoadCalibration )
       {
            return EXIT_FAILURE;
       }
       else  // if calibration information are loaded, go ahead
       {

           // load image filename
           std::vector<std::string> vec_image = stlplus::folder_files( sImageDir );

           // load kept channel
           std::vector< li_Size_t > keptChan;
           loadChannelFile( keptChan, sChannelFile );

           bool isExported = false;

           if( sGpsFileName.empty())
           {
               //create and export list to folder
               isExported =  computeInstrinsicPerImages(
                                       vec_image,
                                       vec_sensorData,
                                       keptChan,
                                       sImageDir,
                                       sOutputDir,
                                       focalPixPermm,
                                       bUseCalibPrincipalPoint,
                                       bRigidRig,
                                       sTimestampLow,
                                       sTimestampUp);
            }
            else
            {
              //create and export list to folder
              isExported =  computeInstrinsicGPSPerImages(
                                      vec_image,
                                      vec_sensorData,
                                      keptChan,
                                      sImageDir,
                                      sOutputDir,
                                      sGpsFileName,
                                      focalPixPermm,
                                      bUseCalibPrincipalPoint,
                                      bRigidRig,
                                      sTimestampLow,
                                      sTimestampUp);

            }

            // do final check to ensure all went well
            if( isExported )
            {
                std::cout << "Sucessfully exported list to folder. Quit" << std::endl;
                return EXIT_SUCCESS;
            }
            else
            {
                std::cerr << "Could not export list to folder. Exit " << std::endl;
                return EXIT_FAILURE;
            }

        }

    }

}
