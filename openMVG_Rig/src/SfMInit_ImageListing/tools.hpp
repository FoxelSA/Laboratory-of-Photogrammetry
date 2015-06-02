/*
 * tools
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

 /*! \file tools.hpp
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

#ifndef TOOLS_HPP_
#define TOOLS_HPP_


#include <fastcal-all.h>
#include <vector>
#include <string.h>
#include <types.hpp>

using namespace std;

/*******************************************************************************
*  Given 4 angles, compute Elphel rotation
*
********************************************************************************
*/

/*! \brief Compute rotation rig referential to sensor referential
*
* This function compute the rotation rig referential to sensor referential using
* elphel calibration angle and rotation.
*
* \param R           Computed rotation
* \param az          Elphel's Angle azimuth (in radian) of subcamera
* \param head        Elphel's Angle heading (in radian) of subcamera
* \param ele         Elphel's Angle elevation (in radian) of subcamera
* \param roll        Elphel's Angle roll (in radian) of subcamera
*
* \return The rotation in the array R
*/


 void computeRotationEl ( li_Real_t* R , li_Real_t az , li_Real_t head, li_Real_t ele , li_Real_t roll);

/********************************************************************************
*  Given three angles, entrance pupil forward, radius and height, compute optical center position.
*
*********************************************************************************
*/

/*! \brief Compute optical center of elphel's subcamera
*
* This function compute the optical center of a given elphel subcamera
*
* \param C           Computed optical center
* \param radius      Radius of optical center in elphel's coordinate frame
* \param heigt       Height of optical center in elphel's coordinate frame
* \param azimuth     Elphel's Angle azimuth (in radian) of subcamera
* \param R           Rotation rig referential frame to sensor frame
* \param entrancePupilForward  Entrance pupil forward of the associated camera
*
* \return The optical center in the array C
*/

 void getOpticalCenter ( li_Real_t* C ,
                const li_Real_t& radius,
                const li_Real_t& height,
                const li_Real_t& azimuth,
                const li_Real_t* R,
                const li_Real_t& entrancePupilForward );

#endif /* TOOLS_HPP_ */
