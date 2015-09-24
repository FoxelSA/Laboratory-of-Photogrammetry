/*
 * relative_translation_test
 *
 * Copyright (c) 2014-2015 FOXEL SA - http://foxel.ch
 * Please read <http://foxel.ch/license> for more information.
 *
 *
 * Author(s):
 *
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

#include "./sfm_robust_relative_pose_rig.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/essential.hpp"
#include "../external/testing/testing.h"

#include "./test_data_sets.hpp"

using namespace openMVG;
using namespace opengv;

/* to check
 * 1) a model is evaluated
 * 2) computed model (difference between rotations + translation )
 * 3) reprojection error is small
 */

TEST(GeRelativePose, GeRelativePose_Kernel) {

  // Use the Generalized eigenvalue solver to solve pose problem
  typedef openMVG::noncentral::kernel::GePointSolver SolverType;

  int focal = 1;
  int principal_Point = 0;

  // initialize structure
  const size_t nPoses = 4;
  const size_t rig_size = 3;
  const size_t nViews = rig_size * nPoses;
  const size_t nbPoints = 2*SolverType::MINIMUM_SAMPLES;

  //-- Setup a circular camera rig and assert that GE relative pose works.
  const NPoseDataSet d = NRealisticPosesRing(nPoses, nViews, rig_size, nbPoints,
    nPoseDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Pose_Estimation.ply");
  //-- Test triangulation of all the point
  NPoseDataSet d2 = d;

  // check that the scene is well constructed
  double  triangulation_error = 0.0;
  for (size_t i = 0; i  < nbPoints; ++i)
  {
    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    std::vector < std::pair <Mat34, Vec2> >  views;

    for(size_t j = 0; j < nPoses; ++j)
      for(size_t k = 0 ; k < rig_size; ++k )
      {
        //compute projection matrices
        Vec2 pt;
        pt << d2._x[j * rig_size + k].col(i)(0), d2._x[j * rig_size + k].col(i)(1);
        views.push_back(std::make_pair ( d.P(j,k), pt ) );
      }

    // update triangulation object
    for( size_t i = 0 ; i < views.size(); ++i )
      triangulationObj.add ( views[i].first, views[i].second );

    const Vec3 X = triangulationObj.compute();
    triangulation_error += (X -d._X.col(i)).norm();
  }

  // Check that triangulation point are near to the inital data
  EXPECT_NEAR(0.0, triangulation_error, 1e-8);

  // initialize opengv offsets and rotations
  opengv::translations_t  offsets;
  opengv::rotations_t     rotations;
  for(int i=0; i < rig_size; ++i )
  {
      offsets.push_back(d._offsets[i]);
      rotations.push_back(d._rotations[i].transpose());
  }

  // evaluate models
  for(int n=0 ; n < nPoses ; ++n )
  {
    // initialize structures used for matching between rigs
    opengv::bearingVectors_t bearingVectorsRigOne, bearingVectorsRigTwo;
    std::vector<int> camCorrespondencesRigOne, camCorrespondencesRigTwo;

    for(int j=0; j < nbPoints; ++j )  // iteration on 3D points
    {
      for(int k=0; k < rig_size; ++k )  // iteration on subcamera of first rig
      {
        for( int l=0; l < rig_size; ++l )  // iterations on subcamera of second rig
        {
          // update sub camera index lists
          camCorrespondencesRigOne.push_back(k);
          camCorrespondencesRigTwo.push_back(l);

          // update bearing vectors lists
          opengv::bearingVector_t  bearing_vector_0;
          opengv::bearingVector_t  bearing_vector_1;

          bearing_vector_0 << d2._x[n*rig_size+k].col(j)(0),
                              d2._x[n*rig_size+k].col(j)(1),
                              1.0;
          bearing_vector_1 << d2._x[(n+1)%nPoses*rig_size+l].col(j)(0),
                              d2._x[(n+1)%nPoses*rig_size+l].col(j)(1),
                              1.0;

          bearing_vector_0.normalize();
          bearing_vector_1.normalize();

          bearingVectorsRigOne.push_back( bearing_vector_0 );
          bearingVectorsRigTwo.push_back( bearing_vector_1 );
        }
      }
    }

    // Define kernel
    typedef openMVG::robust::ACKernelAdaptorRigPose<
           SolverType,
           openMVG::noncentral::kernel::RigAngularError,
           transformation_t>
           KernelType;

    KernelType kernel(  bearingVectorsRigOne,
                        bearingVectorsRigTwo,
                        camCorrespondencesRigOne,
                        camCorrespondencesRigTwo,
                        offsets,
                        rotations);

    std::vector<opengv::transformation_t> models;
    vector<size_t> samples;

    for (size_t k = 0; k < bearingVectorsRigOne.size(); ++k) {
      samples.push_back(k);
    }

    kernel.Fit(samples, &models);
    d2._R[n] = models[0].block<3,3>(0,0).transpose();
    d2._t[n] = -d2._R[n] * models[0].col(3);
    d2._C[n] = models[0].col(3);
  }

  // Assert that found relative motion is correct for almost one model.
  bool bsolution_found = false;
  for (size_t n = 0; n < nPoses; ++n) {
    //-- Compute Ground Truth motion
    Mat3 R;
    Vec3 t;
    RelativeCameraMotion(d._R[n], d._t[n], d._R[(n+1)%nPoses], d._t[(n+1)%nPoses], &R, &t);

    // Check that we find the correct relative orientation.
    Mat3 nul_mat = Mat3::Zero();
    if (FrobeniusDistance(R, d2._R[n]) < 1e-1 * FrobeniusDistance(nul_mat, d2._R[n])
      && (t / t.norm() - d2._t[n] / d2._t[n].norm()).norm() < 1e-1 ) {
        bsolution_found = true;
      }
  }

  //-- Almost one solution must find the correct relative orientation
  CHECK(bsolution_found);
}

TEST(SixPtRelativePoseTest, SixPtRelativePoseTest_Kernel ) {
  // Use the Generalized eigenvalue solver to solve pose problem
  typedef openMVG::noncentral::kernel::SixPointSolver SolverType;

  int focal = 1;
  int principal_Point = 0;

  // initialize structure
  const size_t nPoses = 4;
  const size_t rig_size = 3;
  const size_t nViews = rig_size * nPoses;
  const size_t nbPoints = SolverType::MINIMUM_SAMPLES;

  //-- Setup a circular camera rig and assert that GE relative pose works.
  const NPoseDataSet d = NRealisticPosesRing(nPoses, nViews, rig_size, nbPoints,
  nPoseDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Pose_Estimation.ply");
  //-- Test triangulation of all the point
  NPoseDataSet d2 = d;

  // check that the scene is well constructed
  double  triangulation_error = 0.0;
  for (size_t i = 0; i  < nbPoints; ++i)
  {
    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    std::vector < std::pair <Mat34, Vec2> >  views;

    for(size_t j = 0; j < nPoses; ++j)
      for(size_t k = 0 ; k < rig_size; ++k )
      {
        //compute projection matrices
        Vec2 pt;
        pt << d2._x[j * rig_size + k].col(i)(0), d2._x[j * rig_size + k].col(i)(1);
        views.push_back(std::make_pair ( d.P(j,k), pt ) );
      }

    // update triangulation object
    for( size_t i = 0 ; i < views.size(); ++i )
      triangulationObj.add ( views[i].first, views[i].second );

      const Vec3 X = triangulationObj.compute();
      triangulation_error += (X -d._X.col(i)).norm();
  }

  // Check that triangulation point are near to the inital data
  EXPECT_NEAR(0.0, triangulation_error, 1e-8);

  // initialize opengv offsets and rotations
  opengv::translations_t  offsets;
  opengv::rotations_t     rotations;
  for(int i=0; i < rig_size; ++i )
  {
    offsets.push_back(d._offsets[i]);
    rotations.push_back(d._rotations[i].transpose());
  }

  // evaluate models
  for(int n=0 ; n < nPoses ; ++n )
  {
    // initialize structures used for matching between rigs
    opengv::bearingVectors_t bearingVectorsRigOne, bearingVectorsRigTwo;
    std::vector<int> camCorrespondencesRigOne, camCorrespondencesRigTwo;

    for(int j=0; j < nbPoints; ++j )  // iteration on 3D points
    {
          // update sub camera index lists
          camCorrespondencesRigOne.push_back(1);
          camCorrespondencesRigTwo.push_back(1);

          // update bearing vectors lists
          opengv::bearingVector_t  bearing_vector_0;
          opengv::bearingVector_t  bearing_vector_1;

          bearing_vector_0 << d2._x[n*rig_size+1].col(j)(0),
              d2._x[n*rig_size+1].col(j)(1),
              1.0;
          bearing_vector_1 << d2._x[(n+1)%nPoses*rig_size+1].col(j)(0),
              d2._x[(n+1)%nPoses*rig_size+1].col(j)(1),
              1.0;

          bearing_vector_0.normalize();
          bearing_vector_1.normalize();

          bearingVectorsRigOne.push_back( bearing_vector_0 );
          bearingVectorsRigTwo.push_back( bearing_vector_1 );
    }

    // Define kernel
    typedef openMVG::robust::ACKernelAdaptorRigPose<
        SolverType,
        openMVG::noncentral::kernel::RigAngularError,
        transformation_t>
        KernelType;

    KernelType kernel(  bearingVectorsRigOne,
        bearingVectorsRigTwo,
        camCorrespondencesRigOne,
        camCorrespondencesRigTwo,
        offsets,
        rotations);

    std::vector<opengv::transformation_t> models;
    vector<size_t> samples;

    for (size_t k = 0; k < SolverType::MINIMUM_SAMPLES; ++k) {
      samples.push_back(k);
    }

    kernel.Fit(samples, &models);
    d2._R[n] = models[0].block<3,3>(0,0).transpose();
    d2._t[n] = -d2._R[n] * models[0].col(3);
    d2._C[n] = models[0].col(3);
  }

  // Assert that found relative motion is correct for almost one model.
  bool bsolution_found = false;
  for (size_t n = 0; n < nPoses; ++n) {
    //-- Compute Ground Truth motion
    Mat3 R;
    Vec3 t;
    RelativeCameraMotion(d._R[n], d._t[n], d._R[(n+1)%nPoses], d._t[(n+1)%nPoses], &R, &t);

    // Check that we find the correct relative orientation.
    // Check that we find the correct relative orientation.
    Mat3 nul_mat = Mat3::Zero();
    if (FrobeniusDistance(R, d2._R[n]) < 1e-1 * FrobeniusDistance(nul_mat, d2._R[n])
      && (t / t.norm() - d2._t[n] / d2._t[n].norm()).norm() < 1e-1 ) {
        bsolution_found = true;
      }
  }

  //-- Almost one solution must find the correct relative orientation
//  CHECK(bsolution_found);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
