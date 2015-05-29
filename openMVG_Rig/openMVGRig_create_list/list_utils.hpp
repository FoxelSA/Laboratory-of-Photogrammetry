/*
 * list_utils
 *
 * Copyright (c) 2014-2015 FOXEL SA - http://foxel.ch
 * Please read <http://foxel.ch/license> for more information.
 *
 *
 * Author(s):
 *
 *      Stéphane Flotron <s.flotron@foxel.ch>
 *
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

 /*! \file list_utils.hpp
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


#ifndef LIST_UTILS_HPP_
#define LIST_UTILS_HPP_


#include <fastcal-all.h>
#include <types.hpp>
#include <tools.hpp>
#include <map>
#include <set>
#include <fstream>
#include <sstream>

using namespace std;

/*******************************************************************************
* verify that given timestamp range is valid
*
********************************************************************************
*/

/*! \brief Check input data validity
*
* This function checks that given timestamp range is valid
*
* \param sTimestampLow  Lower bound for timestamp value
* \param sTimestampUp   Upper bound for timestampe value
*
* \return bool value indicating if the timestamp range is valid or not
*/

bool isRangeValid(  const std::string& sTimestampLow,
                    const std::string& sTimestampUp);


/*******************************************************************************
* verify that software argument are valid
*
********************************************************************************
*/

/*! \brief Check input data validity
*
* This function check that the argument given to listomvg are valid.
*
* \param softname       It's argv[0]'s value
* \param sImageDir      The directory containing the image you want to use
* \param sOuputDir      The directory where you want to put your lists.txt file
* \param smacAddress    The mac address of the considered elphel camera
* \param sMountPoint    The mount point of the camera folder
* \param sChannelFile   The file containing the channel you want to use
* \param bRigidRig      Bool value indicating if we use rig structure or not
* \param bUseCalibPrincipalPoint   Bool value indicating if we use (or not) the principal point of calibration
* \param focalPixPermm  the focal length in pixel per mm (optionnal argument)
* \param sTimestampLow  Lower bound for timestamp value
* \param sTimestampUp   Upper bound for timestampe value
*
* \return bool value that says if the arguments are valid or not
*/

bool isInputValid(  const char* softName,
                    const std::string& sImageDir,
                    const std::string& sOutputDir,
                    const std::string& smacAddress,
                    const std::string& sMountPoint,
                    const std::string& sChannelFile,
                    const bool & bRigidRig,
                    const bool & bUseCalibPrincipalPoint,
                    const double & focalPixPermm,
                    const std::string& sTimestampLow,
                    const std::string& sTimestampUp);

/*********************************************************************
 *  load calibration data related to elphel cameras
 *
 *********************************************************************/

 /*! \brief Calibration data loading
 *
 * This function parse channel of considered elephel camera and load all calibration
 * informations needed by openMVG input file.
 *
 * \param vec_sensorData   Vector containing all sensor informations
 * \param sensor_index     the sensor index of elphel camera (between 0 and Channels-1)
 * \param sMountPoint      The mount point of the camera folder
 * \param smacAddress      The mac address of the considered elphel camera
 *
 * \return bool value that says if the loading was sucessfull or not
 */

bool  loadCalibrationData( std::vector< sensorData >  & vec_sensorData,
                    const std::string & sMountPoint,
                    const std::string & smacAddress) ;

/*********************************************************************
 *  read channel file (if exists )
 *
 *********************************************************************/

 /*! \brief Channel File parser
 *
 * This function parse channel file and create a list of suchannel to keep
 *
 * \param keptChan         List of kept channels
 * \param sChannelFile     The complete path of channel file
 *
 * \return The list of kept channel if vector keptChan
 */

void loadChannelFile( std::vector< li_Size_t >  & keptChan,
                    const std::string & sChannelFile );


/*********************************************************************
*  compute camera and rig intrinsic parameters
*
*********************************************************************/

/*! \brief Compute intrinsic needed by openMVG for each image
*
* This function parse image list and create the associated information needed by openMVG
*
* \param vec_image           Image list
* \param vec_sensorData      Calibration data list
* \param keptChan            List of kept channel
* \param sImageDir           Directory containing the images
* \param sOutputDir          Directory containing the lists.txt file
* \param focalPixPermm       Focal length in pixel per mm
* \param bUsePrincipalPoint  bool value indicating if we use (or not) principal point from calibration
* \param bUSeRigidRig        bool value indicating if we use rig structure or not
* \param sTimestampLower     Lower bound for timestamp interval
* \param sTimestampUpper     Upper bound for timestamp interval
*
* \return bool value telling if we have generated the file lists.txt
*/

bool computeInstrinsicPerImages(
                      std::vector<std::string> & vec_image,
                      const std::vector< sensorData > & vec_sensorData,
                      const std::vector< li_Size_t >  & keptChan,
                      const std::string & sImageDir,
                      const std::string & sOutputDir,
                      const double & focalPixPermm,
                      const bool & bUsePrincipalPoint,
                      const bool & bUseRigidRig,
                      std::string& sTimestampLower,
                      std::string& sTimestampUpper);

/*********************************************************************
 *  compute camera and rig intrinsic parameters
 *
 *********************************************************************/

/*! \brief Compute intrinsic needed by openMVG for each image
 *
 * This function parse image list and create the associated information needed by openMVG
 *
 * \param vec_image           Image list
 * \param vec_sensorData      Calibration data list
 * \param keptChan            List of kept channel
 * \param sImageDir           Directory containing the images
 * \param sOutputDir          Directory containing the lists.txt file
 * \param sGpsFile            Filename of the GPS / imu JSON file
 * \param focalPixPermm       Focal length in pixel per mm
 * \param bUsePrincipalPoint  bool value indicating if we use (or not) principal point from calibration
 * \param bUSeRigidRig        bool value indicating if we use rig structure or not
 * \param sTimestampLower     Lower bound for timestamp interval
 * \param sTimestampUpper     Upper bound for timestamp interval
 *
 * \return bool value telling if we have generated the file lists.txt
*/

bool computeInstrinsicGPSPerImages(
            std::vector<std::string> & vec_image,
            const std::vector< sensorData > & vec_sensorData,
            const std::vector< li_Size_t >  & keptChan,
            const std::string & sImageDir,
            const std::string & sOutputDir,
            const std::string & sGpsFile,
            const double & focalPixPermm,
            const bool & bUsePrincipalPoint,
            const bool & bUseRigidRig,
            std::string& sTimestampLower,
            std::string& sTimestampUpper);


/*********************************************************************
 *  compute image intrinsic parameter
 *
 *********************************************************************/

 /*! \brief Compute intrinsic needed by openMVG for a single image
 *
 * This function compute the information needed by openMVG for a single image
 *
 * \param camInfo             Strucutre containing all information needed by openMVG
 * \param vec_sensorData      Calibration data list
 * \param timestamp           Rig timestamp associated to image
 * \param sensor_index        Subchannel number
 * \param focalPixPermm       Focal length in pixel per mm
 * \param bUsePrincipalPoint  bool value indicating if we use (or not) principal point from calibration
 * \param bUSeRigidRig        bool value indicating if we use rig structure or not
 *
 * \return Camera calibration information in the structure camInfo
 */

void computeImageIntrinsic(
                     camInformation & camInfo,
                     const std::vector < sensorData > & vec_sensorData,
                     const std::string & timestamp,
                     const size_t   & sensor_index,
                     const double   & focalPixPermm,
                     const bool     & bUseCalibPrincipalPoint,
                     const bool     & bRigidRig
);


/*********************************************************************
*  keep only most representative rigs
*
*********************************************************************/

/*! \brief Analyze images and keep only most representative rig
*
* This function analyse image list and keep only images related to the
* rig mostly represented. It create a set of image to remove
*
* \param imageToRemove            Set of image to remove
* \param mapSubcamPerTimestamp    Mapping timestamp to associated images and subcameras
* \param imageNumber              The number of input images
* \param sTimestampLower          Lower bound for timestamp interval
* \param sTimestampUpper          Upper bound for timestamp interval
*
* \return The set of image not usable for reconstruction in set imageToRemove
*/

void keepRepresentativeRigs(
           std::set <string> & imageToRemove,
           const std::map<std::string, std::vector<string> > & mapSubcamPerTimestamp,
           const size_t imageNumber,
           const std::string& sTimestampLower,
           const std::string& sTimestampUpper
);

/*********************************************************************
*  write image list and intrinsic to file
*
*********************************************************************/

/*! \brief Export calibration information to file lists.txt
*
* This function exported the calibration information to file lists.txt so
* that openMVG could directly read it for the reconstruction
*
* \param imageToRemove       Set of image to remove
* \param camAndIntrinsics    Set of pair (image name, image intrinsic)
* \param listTXT             The stream that write the information to file lists.txt
* \param bRigidRig           Use rigid rig or not for reconstruction
*
*/

void exportToFile(
          const std::set <string> & imageToRemove,
          const std::set<imageNameAndIntrinsic>& camAndIntrinsics,
          std::ofstream& listTXT,
          const bool & bRigidRig
);


#endif /* LIST_UTILS_HPP_ */
