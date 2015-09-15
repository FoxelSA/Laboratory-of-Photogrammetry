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

#include "./test_data_sets.hpp"
#include "openMVG/numeric/numeric.h"
#include "../external/CppUnitLite/TestHarness.h"
#include "../external/testing/testing.h"

#include "openMVG/multiview/projection.hpp"

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"

#include "openMVG/linearProgramming/bisectionLP.hpp"
#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_From_xi_Ri.hpp"

#include <iostream>
#include <vector>

using namespace openMVG;

using namespace linearProgramming;
using namespace lInfinityCV;

TEST(Translation_Structure_L_Infinity, OSICLP_SOLVER) {

  const size_t nPoses = 3;
  const size_t rig_size = 2;
  const size_t nViews = 6;
  const size_t nbPoints = 12;
  const NPoseDataSet d = NRealisticPosesRing(nPoses, nViews, rig_size, nbPoints,
    nPoseDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity.ply");
  //-- Test triangulation of all the point
  NPoseDataSet d2 = d;

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  Mat megaMat(4, d._n*d._x[0].cols());
  {
    size_t cpt = 0;
    for (size_t i=0; i<d._n;++i)
    {
      const size_t camIndex = i;
      for (size_t j=0; j < d._x[0].cols(); ++j)
      {
        megaMat(0,cpt) = d._x[camIndex].col(j)(0);
        megaMat(1,cpt) = d._x[camIndex].col(j)(1);
        megaMat(2,cpt) = j;
        megaMat(3,cpt) = camIndex;
        cpt++;
      }
    }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution((nViews + nbPoints)*3);

    OSI_CLP_SolverWrapper wrapperOSICLPSolver(vec_solution.size());
    Translation_Structure_L1_ConstraintBuilder cstBuilder( d._R, megaMat);
    EXPECT_TRUE(
      (BisectionLP<Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      wrapperOSICLPSolver,
      cstBuilder,
      &vec_solution,
      1.0,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nViews; ++i)
      {
        size_t index = i*3;
        d2._t[i] = Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i] * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nViews;
        d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i < d2._n; ++i) {
      for(size_t k = 0; k < d._x[0].cols(); ++k)
      {
        xk = Project(d2.P(i), Vec3(d2._X.col(k)));
        xsum += Vec2(( xk - d2._x[i].col(k)).array().pow(2));
      }
    }
    double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4);
  }

  d2.ExportToPLY("test_After_Infinity.ply");
}

TEST(Translation_Structure_L_Infinity, OSICLP_SOLVER_K) {

  const size_t nPoses = 3;
  const size_t rig_size = 2;
  const size_t nViews = 3;
  const size_t nbPoints = 6;
  const NPoseDataSet d = NRealisticPosesRing(nPoses, nViews, rig_size, nbPoints,
    nPoseDatasetConfigurator(1000,1000,500,500,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity.ply");
  //-- Test triangulation of all the point
  NPoseDataSet d2 = d;

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  Mat megaMat(4, d._n*d._x[0].cols());
  {
    size_t cpt = 0;
    for (size_t i=0; i < d._n;++i)
    {
      const size_t camIndex = i;
      for (size_t j=0; j < (size_t)d._x[0].cols(); ++j)
      {
        megaMat(0,cpt) = d._x[camIndex].col(j)(0);
        megaMat(1,cpt) = d._x[camIndex].col(j)(1);
        megaMat(2,cpt) = j;
        megaMat(3,cpt) = camIndex;
        cpt++;
      }
    }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution((nViews + nbPoints)*3);

    std::vector<Mat3> vec_KR(d._R);
    for(int i=0;i < nViews; ++i)
      vec_KR[i] = d._K[0] * d._R[i];

    OSI_CLP_SolverWrapper wrapperOSICLPSolver(vec_solution.size());
    Translation_Structure_L1_ConstraintBuilder cstBuilder( vec_KR, megaMat);
    EXPECT_TRUE(
      (BisectionLP<Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      wrapperOSICLPSolver,
      cstBuilder,
      &vec_solution,
      1.0,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nViews; ++i)
      {
        size_t index = i*3;
        d2._t[i] = d._K[0].inverse() * Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i] * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nViews;
        d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i < d2._n; ++i) {
      for(size_t k = 0; k < (size_t)d._x[0].cols(); ++k)
      {
        xk = Project(d2.P(i), Vec3(d2._X.col(k)));
        xsum += Vec2(( xk - d2._x[i].col(k)).array().pow(2));
      }
    }
    double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4);
  }

  d2.ExportToPLY("test_After_Infinity.ply");
}

#ifdef OPENMVG_HAVE_MOSEK
TEST(Translation_Structure_L_Infinity, MOSEK) {

  const size_t nViews = 3;
  const size_t nbPoints = 6;
  const NPoseDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    NPoseDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity.ply");
  //-- Test triangulation of all the point
  NPoseDataSet d2 = d;

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  Mat megaMat(4, d._n*d._x[0].cols());
  {
    size_t cpt = 0;
    for (size_t i=0; i<d._n;++i)
    {
      const size_t camIndex = i;
      for (size_t j=0; j<d._x[0].cols(); ++j)
      {
        megaMat(0,cpt) = d._x[camIndex].col(j)(0);
        megaMat(1,cpt) = d._x[camIndex].col(j)(1);
        megaMat(2,cpt) = j;
        megaMat(3,cpt) = camIndex;
        cpt++;
      }
    }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution((nViews + nbPoints)*3);

    MOSEK_SolveWrapper wrapperMosek(vec_solution.size());
    Translation_Structure_L1_ConstraintBuilder cstBuilder( d._R, megaMat);
    EXPECT_TRUE(
      (BisectionLP<Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      wrapperMosek,
      cstBuilder,
      &vec_solution,
      1.0,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nViews; ++i)
      {
        size_t index = i*3;
        d2._t[i] = Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i] * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nViews;
        d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i < d2._n; ++i) {
      for(size_t k = 0; k < d._x[0].cols(); ++k)
      {
        xk = Project(d2.P(i), Vec3(d2._X.col(k)));
        xsum += Vec2(( xk - d2._x[i].col(k)).array().pow(2));
      }
    }
    double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4);
  }

  d2.ExportToPLY("test_After_Infinity.ply");
}
#endif // OPENMVG_HAVE_MOSEK

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
