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

#include <types.hpp>
#include <iostream>

using namespace std;

// test equality between intrinsic parameters
  bool operator==(const camInformation& c1, const camInformation& c2)
  {

    bool bSameRot(true), bSameC(true);

    for( size_t i = 0; i < 9 ; ++i )
       bSameRot = ( bSameRot && (c1.R[i] == c2.R[i]) );

    for( size_t j=0 ; j < 3; ++j )
       bSameC = ( bSameC && (c1.C[j] == c2.C[j]) );

    return ( (c1.width == c2.width) && (c1.height == c2.height) && (c1.focal == c2.focal) && (c1.px0 == c2.px0) && (c1.py0 == c2.py0 ) && bSameRot && bSameC );
  }

  bool operator!=(const camInformation& c1, const camInformation& c2)
  {
      return !( c1 == c2);
  }

  bool operator<(const camInformation& c1, const camInformation& c2)
  {
      if( c1.focal < c2.focal )
      {
          return ( c1.focal < c2.focal );
      }
      else
      {
          if( (c1.px0 < c2.px0) )
          {
              return (c1.px0 < c2.px0 );
          }
          else
          {
              return ( c1.C[0] < c2.C[0] );
          }
      }
  }

  bool operator==(const imageNameAndIntrinsic& i1, const imageNameAndIntrinsic& i2) {
      return (i1.first == i2.first);
    }

  bool operator!=(const imageNameAndIntrinsic& i1, const imageNameAndIntrinsic& i2) {
      return !(i1.first == i2.first);
    }

  // Lexicographical ordering of matches. Used to remove duplicates.
  bool operator<(const imageNameAndIntrinsic& i1, const imageNameAndIntrinsic& i2) {
      return i1.first < i2.first;
    }
