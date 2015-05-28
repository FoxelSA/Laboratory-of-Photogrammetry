/*
 * types.hpp
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

 /*! \file types.hpp
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

 #include <fastcal-all.h>
 #include <string>

 using namespace std;

 #ifndef TYPES_HPP_
 #define TYPES_HPP_

 // integers values
    typedef  unsigned int li_Size_t;

 // real values
    typedef  double li_Real_t;

/******************************************************************************
 * sensorData
 *****************************************************************************/

 /*! \struct sensorData
 * \brief structure used to store calibration information
 *
 * This structure is designed to store the needed informations coming from
 * the elphel camera calibration
 *
 * \var sensorData::lfWidth
 *  Width of sensor image
 * \var sensorData::lfHeight
 *  Height of sensor image
 * \var sensorData::lfChannels
 *  Number of channels of elphel camera
 * \var sensorData::lfFocalLength
 *  Focal length in mm
 * \var sensorData::lfPixelSize
 *  pixel size in mm
 * \var sensorData::lfAzimuth
 *  azimuth angle in elphel coordinate frame (in radian)
 * \var sensorData::lfHeading
 *  heading angle in elphel coordinate frame (in radian)
 * \var sensorData::lfElevation
 *  Elevation angle in elphel coordinate frame (in radian)
 * \var sensorData::lfRoll
 *  roll around z axis (in radian)
 * \var sensorData::lfpx0
 *  x coordinate of principal point of sensor image, in pixels
 * \var sensorData::lfpy0
 *  y coordinate of principal point of sensor image, in pixels
 * \var sensorData::lfRadius
 *  radius of optical center of channel in elphel coordinate frame
 * \var sensorData::lfCheight
 *  height of optical center of channel in elphel coordinate frame
 * \var sensorData::lfEntrance
 *  Entrance pupil forward of channel
 * \var sensorData::R
 *  Rotation rig coordinate frame to sensor coordinate frame
 * \var sensorData::C
 *  sensor's optical center in rig coordinate frame
 */

struct sensorData
{
    lf_Size_t   lfWidth     = 0;
    lf_Size_t   lfHeight    = 0;
    lf_Size_t   lfChannels  = 0;

    lf_Real_t   lfFocalLength = 0.0;
    lf_Real_t   lfPixelSize   = 0.0;
    lf_Real_t   lfAzimuth     = 0.0;
    lf_Real_t   lfHeading     = 0.0;
    lf_Real_t   lfElevation   = 0.0;
    lf_Real_t   lfRoll        = 0.0;
    lf_Real_t   lfpx0         = 0.0;
    lf_Real_t   lfpy0         = 0.0;
    lf_Real_t   lfRadius      = 0.0;
    lf_Real_t   lfCheight     = 0.0;
    lf_Real_t   lfEntrance    = 0.0;

    li_Real_t R[9] = {
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0};

    li_Real_t C[3] = {0,0,0};
 };

 /******************************************************************************
 * camera information
 *****************************************************************************/

 /*! \struct camInformation
 * \brief structure used to store calibration information used for file generation
 *
 * This structure is designed to store the needed informations coming from
 * the elphel camera calibration
 *
 * \var camInformation::width
 *  Width of sensor image
 * \var camInformation::height
 *  Height of sensor image
 * \var camInformation::subChan
 *  The subchannel number
 * \var camInformation::focal
 *  Focal length in pixel per mm
 * \var camInformation::px0
 *  X coordinate of principal point
 * \var camInformation::py0
 *  Y coordinate of principal point
 * \var camInformation::R
 *  Rotation rig coordinate frame to sensor coordinate frame
 * \var camInformation::C
 *  sensor's optical center in rig coordinate frame
 */

 struct camInformation
{
  std::string sRigName    = "";

  li_Size_t   width     = 0;
  li_Size_t   height    = 0;
  li_Size_t   subChan   = 0;

  li_Real_t   focal       = 0.0;
  li_Real_t   px0         = 0.0;
  li_Real_t   py0         = 0.0;

  li_Real_t R[9] = {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0};

  li_Real_t C[3] = {0,0,0};
};

/// The structure used to store intrinsic per image
typedef std::pair<std::string, camInformation > imageNameAndIntrinsic;

 // supported image format
 enum Format {
   Pnm, Png, Jpg, Tiff, Unknown
 };

 Format GetFormat(const char *c);

 bool operator==(const imageNameAndIntrinsic& i1, const imageNameAndIntrinsic& i2);

 bool operator!=(const imageNameAndIntrinsic& i1, const imageNameAndIntrinsic& i2);

 // Lexicographical ordering of matches. Used to remove duplicates.
 bool operator<(const imageNameAndIntrinsic& i1, const imageNameAndIntrinsic& i2);

 #endif /* TYPES_HPP_ */
