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

#include <cmath>
#include <iostream>
#include <fstream>

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/projection.hpp"
#include "./test_data_sets.hpp"

namespace openMVG {


nPoseDatasetConfigurator::nPoseDatasetConfigurator(int fx, int fy,
  int cx, int cy, double distance, double jitter_amount):
  _fx(fx), _fy(fy), _cx(cx), _cy(cy), _dist(distance),
  _jitter_amount(jitter_amount)
{}

NPoseDataSet NRealisticPosesRing(  size_t nposes, size_t nviews,
                                   size_t rig_size, size_t npoints,
                                   const nPoseDatasetConfigurator config)
{
  //-- Setup a camera circle rig.
  NPoseDataSet d;
  d._n = nposes;
  d._K.resize(rig_size);
  d._offsets.resize(rig_size);
  d._rotations.resize(rig_size);
  d._R.resize(nposes);
  d._t.resize(nposes);
  d._C.resize(nposes);
  d._x.resize(nviews);
  d._x_ids.resize(nviews);

  d._X.resize(3, npoints);
  d._X.setRandom();
  d._X *= 0.6;

  Vecu all_point_ids(npoints);
  // initialize point id
  for (size_t j = 0; j < npoints; ++j)
    all_point_ids[j] = j;

  // initialize camera rig intrinsic parameters
  for (size_t j = 0 ; j < rig_size ; ++j )
  {
      Vec3  center, lookdir;

      // initialize camera matrix
      d._K[j] << config._fx,           0, config._cx,
                        0,  config._fy, config._cy,
                        0,           0,          1;

      // intialize offset and rotations
      double theta = j * 2 * M_PI / rig_size / nposes ;
      center  << -sin(theta), 0.0, -cos(theta);
      center *= config._dist;
      if( j < 1)
        d._offsets[j] = center ;
      else
        d._offsets[j] = center - d._offsets[0];

      if( j == rig_size - 1 )
        d._offsets[0] *= 0.0;

      // initialize sub camera rotation
      lookdir = -center;
      d._rotations[j] = LookAt(lookdir);
  }

  // now initialize 2d correspondences
  for (size_t i = 0; i < nposes; ++i) {
    Vec3 pose_center, t, jitter, lookdir;

    const double theta = i * 2 * M_PI / nposes;
    //-- Circle equation
    pose_center << sin(theta), 0.0, cos(theta); // Y axis UP
    pose_center *= config._dist;
    d._C[i] = pose_center;

    jitter.setRandom();
    jitter *= config._jitter_amount / pose_center.norm();
    lookdir = -pose_center + jitter;

    d._R[i] = LookAt(lookdir);  // Y axis UP
    d._t[i] = -d._R[i] * pose_center; // [t]=[-RC] Cf HZ.
    for( size_t j = 0; j < rig_size; ++j )
    {
      d._x[i * rig_size + j] = Project(d.P(i, j), d._X);
      d._x_ids[i * rig_size + j] = all_point_ids;
    }
  }
  return d;
}

Mat34 NPoseDataSet::P(size_t i, size_t j)const {
  assert(i < _n);
  assert(j < _K.size());
  Mat34 P;
  P_From_KRt(_K[j], _rotations[j]*_R[i], _rotations[j] * _t[i] - _rotations[j] * _offsets[j], &P);
  return P;
}

void NPoseDataSet::ExportToPLY(
  const std::string & out_file_name)const {
  std::ofstream outfile;
  outfile.open(out_file_name.c_str(), std::ios_base::out);
  if (outfile.is_open()) {
    outfile << "ply"
     << std::endl << "format ascii 1.0"
     << std::endl << "comment NPoseDataSet export"
     << std::endl << "comment It shows 3D point structure and cameras"
                  << "+ camera looking direction"
     << std::endl << "element vertex " << _X.cols() + _t.size()*_K.size()*2
     << std::endl << "property float x"
     << std::endl << "property float y"
     << std::endl << "property float z"
     << std::endl << "property uchar red"
     << std::endl << "property uchar green"
     << std::endl << "property uchar blue"
     << std::endl << "end_header" << std::endl;

    //-- Export 3D point cloud
    for(Mat3X::Index i = 0; i < _X.cols(); ++i) {
      // Exports the point position and point color
      outfile << _X.col(i).transpose()
        << " " << "255 255 255" << std::endl;
    }

    //-- Export 3D camera position t = -RC
    for(size_t i = 0; i < _t.size(); ++i) {
      for(size_t j = 0; j < _K.size(); ++j)
      {
        const Vec3 center = -_R[i].transpose() * _rotations[j].transpose() * ( _rotations[j] * _t[i] - _rotations[j] * _offsets[j] );
        // Exports the camera position and camera color
        outfile << center.transpose() << " " << "0 255 0" << std::endl;
      }
    }
    for(size_t i = 0; i < _t.size(); ++i) {
      for(size_t j = 0; j < _K.size(); ++j)
      {
        Vec3 test;
        test << 0, 0 , 0.4;

        const Vec3 center = -_R[i].transpose() * _rotations[j].transpose() * ( _rotations[j] * _t[i] - _rotations[j] * _offsets[j] );
        // Exports the camera normal
        outfile << center.transpose() +
          (_R[i].transpose()*_rotations[j].transpose()*test).transpose()
          << " " << "255 0 0" << std::endl;
      }
    }
    outfile.close();
  }
}

}  // namespace openMVG
