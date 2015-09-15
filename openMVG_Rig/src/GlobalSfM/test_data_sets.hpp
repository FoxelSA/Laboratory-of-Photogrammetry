/*
 * relative_translation_test
 *
 * Copyright (c) 2014-2015 FOXEL SA - http://foxel.ch
 * Please read <http://foxel.ch/license> for more information.
 *
 *
 * Author(s):
 *
 *      Pierre Moulon <p.moulon@foxel.ch>
 *      Stephane Flotron <s.flotron@foxel.ch>
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

#ifndef OPENMVG_MULTIVIEW_TEST_DATA_SETS_H_
#define OPENMVG_MULTIVIEW_TEST_DATA_SETS_H_

#include "openMVG/numeric/numeric.h"
#include <vector>

namespace openMVG {

using namespace std;

// A N-view metric dataset.
// All points are seen by all cameras.
struct NViewDataSet {
  vector<Mat3> _K;   // Internal parameters (fx, fy, etc).
  vector<Mat3> _R;   // Rotation.
  vector<Vec3> _t;   // Translation.
  vector<Vec3> _C;   // Camera centers.
  Mat3X _X;          // 3D points.
  vector<Mat2X> _x;  // Projected points; may have noise added.
  vector<Vecu>  _x_ids;// Indexes of points corresponding to the projections

  size_t _n;  // Actual number of cameras.

  //-- Return P=K*[R|t] for the Inth camera
  Mat34 P(size_t i) const;

  /// Export in PLY the point structure and camera and camera looking dir.
  void ExportToPLY(const std::string & out_file_name) const;
};

struct nViewDatasetConfigurator
{
  /// Internal camera parameters (focal, principal point)
  int _fx, _fy, _cx, _cy;

  /// Camera random position parameters
  double _dist;
  double _jitter_amount;

  nViewDatasetConfigurator(int fx = 1000,  int fy = 1000,
                           int cx = 500,   int cy  = 500,
                           double distance = 1.5,
                           double jitter_amount = 0.01 );
};

/// Place cameras on a circle with point in the center
NViewDataSet NRealisticCamerasRing(size_t nviews, size_t npoints,
                                   const nViewDatasetConfigurator
                                     config = nViewDatasetConfigurator());

/// Place cameras on cardiod shape with point in the center
NViewDataSet NRealisticCamerasCardioid(size_t nviews, size_t npoints,
                                       const nViewDatasetConfigurator
                                        config = nViewDatasetConfigurator());

} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TEST_DATA_SETS_H_
