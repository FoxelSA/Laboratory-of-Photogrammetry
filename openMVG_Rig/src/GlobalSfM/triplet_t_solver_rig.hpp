/*
 * triplet_t_solver
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

#pragma once

#include "openMVG/numeric/numeric.h"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include <fstream>
#include <utility>
#include <vector>

namespace openMVG   {
namespace lInfinityCV  {

//-- Estimate the translation and the structure
//    from known image points coordinates and camera rotations.
//
// Compare to OpenMVG implementation, the following implementation
//  presents an extension to rig of cameras (camera with a known subposes).
//

/// Encode translation and structure linear program for camera with known subposes
void EncodeRigTiXi
(
  const Mat & M, //Scene representation
  const std::vector<Mat3> Ri,
  const std::vector<Mat3> rigRotation,
  const std::vector<Vec3> rigOffsets,
  double sigma, // Start upper bound
  sRMat & A,
  Vec & C,
  std::vector<openMVG::linearProgramming::LP_Constraints::eLP_SIGN> & vec_sign,
  std::vector<double> & vec_costs,
  std::vector< std::pair<double,double> > & vec_bounds
)
{
  // Build Constraint matrix.
  const size_t Nrig = (size_t) M.row(4).maxCoeff()+1;
  const size_t N3D  = (size_t) M.row(2).maxCoeff()+1;
  const size_t Nobs = M.cols();

  assert(Nrig == Ri.size());

  A.resize(5 * Nobs, 3 * (N3D + Nrig) );

  C.resize(5 * Nobs, 1);
  C.fill(0.0);
  vec_sign.resize(5 * Nobs + 3);

  const size_t transStart  = 0;
  const size_t pointStart  = transStart + 3*Nrig;

# define TVAR(i, el) (0 + 3*(i) + (el))
# define XVAR(j, el) (pointStart + 3*(j) + (el))

  // By default set free variable:
  vec_bounds = std::vector< std::pair<double,double> >(3 * (N3D + Nrig),
    std::make_pair((double)-1e+30, (double)1e+30));
  // Fix the translation ambiguity. (set first cam at (0,0,0))
  vec_bounds[0] = vec_bounds[1] = vec_bounds[2] = std::make_pair(0,0);

  size_t rowPos = 0;
  // Add the cheirality conditions (R_c*R_i*X_j + R_c*T_i + t_c)_3 + Z_ij >= 1
  for (size_t k = 0; k < Nobs; ++k)
  {
    const size_t indexPt3D = M(2,k);
    const size_t indexCam  = M(3,k);
    const size_t indexRig  = M(4,k);

    const Mat3 & R  = Ri[indexRig];
    const Mat3 & Rc = rigRotation[indexCam];
    const Vec3 tc = -Rc * rigOffsets[indexCam];

    const Mat3 RcRi = Rc * R;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(2,2);
    C(rowPos) = 1e-2 - tc(2); // Force a minimum depth to be at least 1e-2 meters
    vec_sign[rowPos] = openMVG::linearProgramming::LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    const Vec2 pt  = M.block<2,1>(0,k);
    const double u = pt(0);
    const double v = pt(1);

    // x-residual =>
    // (R_c*R_i*X_j + R_c*T_i + T_c)_1 / (R_c*R_i*X_j + R_c*T_i + T_c)_3 - u >= -sigma
    // (R_c*R_i*X_j + R_c*T_i + T_c)_1 - u * (R_c*R_i*X_j + R_c*T_i + T_c)_3  + sigma (R_c*R_i*X_j + R_c*T_i + T_c)_3  >= 0.0
    // ((R_c*R_i)_3 * (sigma-u) + (R_c*R_i)_1) * X_j +
    //     + (R_c_3 * (sigma-u) + R_c_1)*t_i + (t_c_1 + t_c_3 * (sigma-u) ) >= 0

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(0,0) + (sigma-u) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(0,1) + (sigma-u) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(0,2) + (sigma-u) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(0,0) + (sigma-u) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(0,1) + (sigma-u) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(0,2) + (sigma-u) * Rc(2,2);
    C(rowPos) = -tc(0) -tc(2) * (sigma-u);
    vec_sign[rowPos] = openMVG::linearProgramming::LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(0,0) - (sigma+u) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(0,1) - (sigma+u) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(0,2) - (sigma+u) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(0,0) - (sigma+u) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(0,1) - (sigma+u) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(0,2) - (sigma+u) * Rc(2,2);
    C(rowPos) = -tc(0) + tc(2) * (sigma + u);
    vec_sign[rowPos] = openMVG::linearProgramming::LP_Constraints::LP_LESS_OR_EQUAL;
    ++rowPos;

    // y-residual =>
    // (R_c*R_i*X_j + R_c*T_i + T_c)_2 / (R_c*R_i*X_j + R_c*T_i + T_c)_3 - v >= -sigma
    // (R_c*R_i*X_j + R_c*T_i + T_c)_2 - v * (R_c*R_i*X_j + R_c*T_i + T_c)_3  + sigma (R_c*R_i*X_j + R_c*T_i + T_c)_3  >= 0.0
    // ((R_c*R_i)_3 * (sigma-v) + (R_c*R_i)_2) * X_j +
    //     + (R_c_3 * (sigma-v) + R_c_2)*t_i + (t_c_2 + t_c_3 * (sigma-v) ) >= 0

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(1,0) + (sigma-v) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(1,1) + (sigma-v) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(1,2) + (sigma-v) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(1,0) + (sigma-v) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(1,1) + (sigma-v) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(1,2) + (sigma-v) * Rc(2,2);
    C(rowPos) = -tc(1) -tc(2) * (sigma-v);
    vec_sign[rowPos] = openMVG::linearProgramming::LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(1,0) - (sigma+v) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(1,1) - (sigma+v) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(1,2) - (sigma+v) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(1,0) - (sigma+v) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(1,1) - (sigma+v) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(1,2) - (sigma+v) * Rc(2,2);
    C(rowPos) = -tc(1) + tc(2) * (sigma+v);
    vec_sign[rowPos] = openMVG::linearProgramming::LP_Constraints::LP_LESS_OR_EQUAL;
    ++rowPos;
  }
# undef TVAR
# undef XVAR
}

/// Encode translation and structure linear program
void EncodeRigCiXi
(
  const Mat & M, //Scene representation
  const std::vector<Mat3> Ri,
  const std::vector<Mat3> rigRotation,
  const std::vector<Vec3> rigOffsets,
  double sigma, // Start upper bound
  sRMat & A,
  Vec & C,
  std::vector<openMVG::linearProgramming::LP_Constraints::eLP_SIGN> & vec_sign,
  std::vector<double> & vec_costs,
  std::vector< std::pair<double,double> > & vec_bounds
)
{
  // Build Constraint matrix.
  const size_t Nrig = (size_t) M.row(4).maxCoeff()+1;
  const size_t N3D  = (size_t) M.row(2).maxCoeff()+1;
  const size_t Nobs = M.cols();

  assert(Nrig == Ri.size());

  A.resize( 3 * Nobs * (Nobs-1), Nobs + 3*Nrig );

  C.resize( 3 * Nobs * (Nobs-1), 1);
  C.fill(0.0);
  vec_sign.resize(3 * Nobs * (Nobs-1));

  const size_t transStart = 0;
  const size_t pointStart = transStart + 3*Nrig;

# define TVAR(i, el) (0 + 3*(i) + (el))           // translation between rigs
# define XVAR(j) (pointStart + (j))

  // By default set free variable:
  vec_bounds = std::vector< std::pair<double,double> >(Nobs + 3*Nrig,
    std::make_pair((double)-1e+30, (double)1e+30));
  // Fix the translation ambiguity. (set first cam at (0,0,0))
  vec_bounds[0] = vec_bounds[1] = vec_bounds[2] = std::make_pair(0,0);

  Vec3 b0, b1;
  // Add the cheirality conditions (R_c*R_i*X_j + R_c*T_i + t_c)_3 + Z_ij >= 1
  for (size_t k = 0; k < Nobs ; ++k)
  {
    vec_bounds[XVAR(k)] = std::make_pair( 0.0 , (double)1e+30 );
  }

  size_t rowpos = 0;
  // Add the cheirality conditions (R_c*R_i*X_j + R_c*T_i + t_c)_3 + Z_ij >= 1
  for (size_t k = 0; k < Nobs-1 ; ++k)
  {
    for (size_t l = k+1; l < Nobs; ++l)
    {
      // define pose index
      const size_t  pose_I = k;
      const size_t  pose_J = l;

      // we assume here that each track is of length 3
      // extract bearing vectors
      b0 << M(0, pose_I), M(1, pose_I), 1.0;
      b1 << M(0, pose_J), M(1, pose_J), 1.0;

      // extract sub poses rotations (rotation sensor to pose referential )
      const Mat3  Rc0  = rigRotation[M(3, pose_I)].transpose();
      const Mat3  Rc1  = rigRotation[M(3, pose_J)].transpose();

      // extract rotations of poses (rotation poses referential to world referential )
      const Mat3  R0  = Ri[M(4, pose_I)].transpose();
      const Mat3  R1  = Ri[M(4, pose_J)].transpose();

      // compute bearing vector in world referential frame
      // just apply rotation of subpose-> pose and pose -> world
      const Vec3  Rb0 = R0 * Rc0 * b0;
      const Vec3  Rb1 = R1 * Rc1 * b1;

      // compute center of camera in world referential frame
      // apply rotation pose -> world at each offset
      const Vec3  R_c0 = R0 * rigOffsets[M(3, pose_I)];
      const Vec3  R_c1 = R1 * rigOffsets[M(3, pose_J)];

      // 3D point index
      const size_t  pointIndex_I = M(2, pose_I);
      const size_t  pointIndex_J = M(2, pose_J);

      // pose index
      const size_t  pose_index_I = M(4, pose_I);
      const size_t  pose_index_J = M(4, pose_J);

      /**************************************************************
      * a 3d point originated from a bearing vector \b
      * of a subcamera with pose (R_c, C_c) in a rig with pose (R_r, C_r)
      *  is given by the following equation
      *
      *    X = \alpha R_r^T  R_c^T \b - R_r^T t_r + R_r^T C_c
      *
      * where \alpha is the depth of point X related to subcamera and
      * \t_r is the translation of the rig. In this scheme we compute
      * three 3D points X_0, X_1 and X_2 originated from poses 0, 1, 2
      * and translation \t_0, \t_1 and \t_2 such that
      *
      *   \| X_0 - X_1 \|_{\infty}  \leq \sigma
      *
      ****************************************************************
      */
      if( pointIndex_J == pointIndex_I && pose_index_I != pose_index_J )
      {
        // encode matrix
        for( int i=0 ; i < 3 ; ++i )  // loop on componant of translation
        {
          // ||X_0 -X_1 || \leq \sigma is equivalent to
          //  \alpha_0 (R_r0 ^T R_c0^T \b0)_i -  \alpha_1 (R_r1 ^T R_c1^T \b1)_i
          //       -(R_r0^T \t_0)_i + (R_r1^T \t_1)_i
          //        \leq \sigma - R_r0^T C_c0 + R_r1^T C_c1
          A.coeffRef(rowpos, TVAR(pose_index_I, 0)) = -R0(i,0);
          A.coeffRef(rowpos, TVAR(pose_index_I, 1)) = -R0(i,1);
          A.coeffRef(rowpos, TVAR(pose_index_I, 2)) = -R0(i,2);
          A.coeffRef(rowpos, TVAR(pose_index_J, 0)) =  R1(i,0);
          A.coeffRef(rowpos, TVAR(pose_index_J, 1)) =  R1(i,1);
          A.coeffRef(rowpos, TVAR(pose_index_J, 2)) =  R1(i,2);
          A.coeffRef(rowpos, XVAR(pose_I)) =  Rb0(i);
          A.coeffRef(rowpos, XVAR(pose_J)) = -Rb1(i);
          C(rowpos) = sigma - R_c0(i) + R_c1(i);
          vec_sign[rowpos] = openMVG::linearProgramming::LP_Constraints::LP_LESS_OR_EQUAL;
          ++rowpos;

          // ||X_0 -X_1 || \leq \sigma is equivalent to
          //  \alpha_0 (R_r0 ^T R_c0^T \b0)_i -  \alpha_1 (R_r1 ^T R_c1^T \b1)_i
          //       -(R_r0^T \t_0)_i + (R_r1^T \t_1)_i
          //        \geq - \sigma - R_r0^T C_c0 + R_r1^T C_c1
          A.coeffRef(rowpos, TVAR(pose_index_I, 0)) = -R0(i,0);
          A.coeffRef(rowpos, TVAR(pose_index_I, 1)) = -R0(i,1);
          A.coeffRef(rowpos, TVAR(pose_index_I, 2)) = -R0(i,2);
          A.coeffRef(rowpos, TVAR(pose_index_J, 0)) =  R1(i,0);
          A.coeffRef(rowpos, TVAR(pose_index_J, 1)) =  R1(i,1);
          A.coeffRef(rowpos, TVAR(pose_index_J, 2)) =  R1(i,2);
          A.coeffRef(rowpos, XVAR(pose_I)) =  Rb0(i);
          A.coeffRef(rowpos, XVAR(pose_J)) = -Rb1(i);
          C(rowpos) = -sigma - R_c0(i) + R_c1(i);
          vec_sign[rowpos] = openMVG::linearProgramming::LP_Constraints::LP_GREATER_OR_EQUAL;
          ++rowpos;
        }
      }
    }
  }
# undef TVAR
# undef XVAR
}

/// Kernel that set Linear constraints for the
///   - Translation Registration and Structure Problem.
///  Designed to be used with bisectionLP and LP_Solver interface.
///
/// Solve the "Translation Registration and Structure Problem"
///  for 'rig' cameras with known rotations by using a sparse Linear Program.
/// Note: Rig camera: camera with known subposes.

struct Rig_Translation_Structure_L1_ConstraintBuilder
{
  Rig_Translation_Structure_L1_ConstraintBuilder(
    const std::vector<Mat3> & vec_Ri,
    const Mat & M,
    const std::vector<Mat3> & rigRotation,
    const std::vector<Vec3> & rigOffsets):
    _M(M),
    _vec_Ri(vec_Ri),
    _rigRotation(rigRotation),
    _rigOffsets(rigOffsets)
  {
  }

  /// Setup constraints for the translation and structure problem,
  ///  in the openMVG::linearProgramming::LP_Constraints object.
  bool Build(double gamma, openMVG::linearProgramming::LP_Constraints_Sparse & constraint)
  {
    EncodeRigTiXi(
      _M,
      _vec_Ri,
      _rigRotation,
      _rigOffsets,
      gamma,
      constraint._constraintMat,
      constraint._Cst_objective,
      constraint._vec_sign,
      constraint._vec_cost,
      constraint._vec_bounds);

    //-- Setup additional information about the Linear Program constraint
    // We look for nb translations and nb 3D points.
    const size_t N3D  = (size_t) _M.row(2).maxCoeff() + 1;
    const size_t Nrig = (size_t) _M.row(4).maxCoeff() + 1;

    constraint._nbParams = (Nrig + N3D) * 3;

    return true;
  }

  std::vector<Mat3> _vec_Ri;  // Rotation matrix
  Mat _M; // M contains (X,Y,index3dPoint, indexCam)^T
  std::vector<Mat3> _rigRotation; // rotation of rig subcameras
  std::vector<Vec3> _rigOffsets; // optical center of rig subcameras in rig referential frame
};

/// Kernel that set Linear constraints for the
///   - Translation Registration and Structure Problem.
///  Designed to be used with bisectionLP and LP_Solver interface.
///
/// Solve the "Translation Registration and Structure Problem"
///  for 'rig' cameras with known rotations by using a sparse Linear Program.
/// Note: Rig camera: camera with known subposes.

struct Rig_Center_Structure_L1_ConstraintBuilder
{
  Rig_Center_Structure_L1_ConstraintBuilder(
    const std::vector<Mat3> & vec_Ri,
    const Mat & M,
    const std::vector<Mat3> & rigRotation,
    const std::vector<Vec3> & rigOffsets):
    _M(M),
    _vec_Ri(vec_Ri),
    _rigRotation(rigRotation),
    _rigOffsets(rigOffsets)
  {
  }

  /// Setup constraints for the translation and structure problem,
  ///  in the openMVG::linearProgramming::LP_Constraints object.
  bool Build(double gamma, openMVG::linearProgramming::LP_Constraints_Sparse & constraint)
  {
    EncodeRigCiXi(
      _M,
      _vec_Ri,
      _rigRotation,
      _rigOffsets,
      gamma,
      constraint._constraintMat,
      constraint._Cst_objective,
      constraint._vec_sign,
      constraint._vec_cost,
      constraint._vec_bounds);

    //-- Setup additional information about the Linear Program constraint
    // We look for nb translations and nb 3D points.
    const size_t N3D  = (size_t) _M.row(2).maxCoeff() + 1;
    const size_t Nobs = (size_t) _M.cols();
    const size_t Nrig = (size_t) _M.row(4).maxCoeff() + 1;

    constraint._nbParams = Nobs + Nrig* 3;

    return true;
  }

  std::vector<Mat3> _vec_Ri;  // Rotation matrix
  Mat _M; // M contains (X,Y,index3dPoint, indexCam)^T
  std::vector<Mat3> _rigRotation; // rotation of rig subcameras
  std::vector<Vec3> _rigOffsets; // optical center of rig subcameras in rig referential frame
};

} // namespace lInfinityCV
} // namespace openMVG
