/*
 * sfm_robust_relative_pose_rig
 *
 * Copyright (c) 2014-2015 FOXEL SA - http://foxel.ch
 * Please read <http://foxel.ch/license> for more information.
 *
 *
 * Author(s):
 *
 *      Stephane Flotron  <s.flotron@foxel.ch>
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

#include <vector>
#include "openMVG/numeric/numeric.h"

#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/cameras/PinholeCamera.hpp"

#include "opengv/types.hpp"
#include "opengv/relative_pose/methods.hpp"
#include "opengv/triangulation/methods.hpp"
#include "opengv/relative_pose/NoncentralRelativeAdapter.hpp"
#include "opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp"

/*************************
 *
 *  Definition of AC ransac kernel adaptor for rig relative poses
 *
 ************************
  */

namespace openMVG {
  namespace robust{

   using namespace opengv;
    /// Rig Pose Kernel adaptator for the A contrario model estimator
    template <typename SolverArg,
      typename ErrorArg,
      typename ModelArg = transformation_t >
    class ACKernelAdaptorRigPose
    {
    public:
      typedef SolverArg Solver;
      typedef ModelArg  Model;
      typedef ErrorArg ErrorT;

      ACKernelAdaptorRigPose(
        const bearingVectors_t & b1,
        const bearingVectors_t & b2,
        const std::vector<int> & scIdOne,
        const std::vector<int> & scIdTwo,
        const translations_t & rigOffsets,
        const rotations_t & rigRotations )
        : b1_(b1), b2_(b2),
          scIdOne_(scIdOne), scIdTwo_(scIdTwo),
          offSets(rigOffsets), rotations(rigRotations),
          logalpha0_(log10(M_PI))
    {
        assert(b1_.size() == b2_.size());
        logalpha0_ = log10(0.5);
      }

      enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
      enum { MAX_MODELS = Solver::MAX_MODELS };

      void Fit(const std::vector<size_t> &samples, std::vector<Model> *models)
      const {
        //create non-central relative adapter
        relative_pose::NoncentralRelativeAdapter adapter(
              b1_, b2_, scIdOne_ , scIdTwo_,
              offSets, rotations);

        Solver::Solve(adapter, models, samples);
      }

      void Errors(const Model & model, std::vector<double> & vec_errors) const
      {
        //create non-central relative adapter
        relative_pose::NoncentralRelativeAdapter adapter(
              b1_, b2_, scIdOne_ , scIdTwo_,
              offSets, rotations);

        vec_errors.resize(b1_.size());
        for (size_t sample = 0; sample < b1_.size(); ++sample)
          vec_errors[sample] = ErrorT::Error(sample, model, adapter);
      }

      double Error(size_t sample, const Model & model) const
      {
        //create non-central relative adapter
        relative_pose::NoncentralRelativeAdapter adapter(
              b1_, b2_, scIdOne_ , scIdTwo_,
              offSets, rotations);

        return ErrorT::Error(sample, model, adapter);
      }

      size_t NumSamples() const { return b1_.size(); }
      void Unnormalize(Model * model) const {}
      double logalpha0() const {return logalpha0_;}
      double multError() const {return 0.5;} // angle vs. angle error
      Mat3 normalizer1() const {return Mat3::Identity();}
      Mat3 normalizer2() const {return Mat3::Identity();}
      double unormalizeError(double val) const { return val; }

    private:
      bearingVectors_t b1_, b2_; // bearing vectors.
      std::vector<int> scIdOne_, scIdTwo_ ; // cam correspondences
      translations_t offSets ; // camera position in rigid rig frame
      rotations_t  rotations ; // rig subcamera orientation
      double logalpha0_; // Alpha0 is used to make the error adaptive to the image size
    };
  } // end ROBUST
} // end openMVG

namespace openMVG {
namespace noncentral {
namespace kernel {

using namespace std;
using namespace opengv;
using namespace openMVG;

/**
 * Six point solver for non central camera system,
 * // [1] "Solutions to minimal generalized relative pose problems".
 * // authors: Stewenius, H., Nister, D., Oskarsson, M., & Astrom K,
 * // Date: 2005:
 * // Conference: Workshop on omnidirectional vision 2005.
 */
struct SixPointSolver {
  enum { MINIMUM_SAMPLES = 9 };
  enum { MAX_MODELS = 64 };
  static void Solve(relative_pose::NoncentralRelativeAdapter & adapter,
                    std::vector<transformation_t> * models,
                    const std::vector<size_t> &indices);
};

/**
 * Generalized eigen value solver for non central camera system,
 * @InProceedings{ Kneip_2014_CVPR,
 * author = {Kneip, Laurent and Li, Hongdong},
 * title = {Efficient Computation of Relative Pose for Multi-Camera Systems},
 * journal = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
 * month = {June},
 * year = {2014}
 * }
 */

struct GePointSolver {
  enum { MINIMUM_SAMPLES = 6 };
  enum { MAX_MODELS = 1 };
  static void Solve(relative_pose::NoncentralRelativeAdapter & adapter,
  std::vector<transformation_t> * models,
  const std::vector<size_t> &indices);
};

// compute reprojection error
struct RigProjError {
  static double Error(size_t sample,
                      const transformation_t & relativePose,
                      relative_pose::NoncentralRelativeAdapter & _adapter)
  {
    // extract pose of cameras
    const Vec3 bearingOne = _adapter.getBearingVector1(sample);
    const Vec3 bearingTwo = _adapter.getBearingVector2(sample);

    const Vec2 x1 = bearingOne.head(2) / bearingOne(2);
    const Vec2 x2 = bearingTwo.head(2) / bearingTwo(2);

    const Mat3 R1 = _adapter.getCamRotation1(sample).transpose();
    const Mat3 R2 = _adapter.getCamRotation2(sample).transpose();

    const Vec3 t1 = - R1 * _adapter.getCamOffset1(sample);
    const Vec3 t2 = - R2 * _adapter.getCamOffset2(sample);

    // retrieve relative pose of rigs
    const translation_t CRig = relativePose.col(3);
    const rotation_t rRig = relativePose.block<3,3>(0,0).transpose();
    const Vec3  tRig = -rRig * CRig;

    // compute relative pose of cameras
    const rotation_t R = R2 * rRig ;
    const translation_t t = R2 * tRig + t2 ;

    // compute 3d point and reprojection error
    const Mat34 P1 = HStack(R1, t1);
    const Mat34 P2 = HStack(R, t);

    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    triangulationObj.add(P1, x1);
    triangulationObj.add(P2, x2);

    const Vec3 X = triangulationObj.compute();

    //- Return max error as a test
    double pt1ReProj = (Project(P1, X) - x1).norm();
    double pt2ReProj = (Project(P2, X) - x2).norm();

    return std::max(pt1ReProj, pt2ReProj);
  }
};
typedef RigProjError SimpleError;

// compute angular error (as in kneip opengv library)
struct RigAngularError {
  static double Error(size_t index,
                      const transformation_t & model,
                      relative_pose::NoncentralRelativeAdapter & _adapter)
  {
    // extract rotation and translation from model
    const translation_t translation = model.col(3);
    const rotation_t rotation = model.block<3,3>(0,0);

    // compute pose
    const translation_t cam1Offset = _adapter.getCamOffset1(index);
    const rotation_t cam1Rotation = _adapter.getCamRotation1(index);
    const translation_t cam2Offset = _adapter.getCamOffset2(index);
    const rotation_t cam2Rotation = _adapter.getCamRotation2(index);

    const translation_t directTranslation =
        cam1Rotation.transpose() *
        ((translation - cam1Offset) + rotation * cam2Offset);
    const rotation_t directRotation =
        cam1Rotation.transpose() * rotation * cam2Rotation;

    _adapter.sett12(directTranslation);
    _adapter.setR12(directRotation);

    transformation_t inverseSolution;
    inverseSolution.block<3,3>(0,0) = directRotation.transpose();
    inverseSolution.col(3) =
        -inverseSolution.block<3,3>(0,0)*directTranslation;

    Vec4 p_hom;
    p_hom << opengv::triangulation::triangulate2(_adapter,index), 1.0;

    const bearingVector_t reprojection1 = p_hom.block<3,1>(0,0).normalized();
    const bearingVector_t reprojection2 = (inverseSolution * p_hom).normalized();
    const bearingVector_t f1 = _adapter.getBearingVector1(index).normalized();
    const bearingVector_t f2 = _adapter.getBearingVector2(index).normalized();

    //bearing-vector based outlier criterium (select threshold accordingly):
    //1-(f1'*f2) = 1-cos(alpha) \in [0:2]
    const double reprojError1 = 1.0 - f1.transpose() * reprojection1 ;
    const double reprojError2 = 1.0 - f2.transpose() * reprojection2 ;
    return std::max(reprojError1,reprojError2);
  }
};
typedef RigAngularError SimpleAngularError;

}  // namespace kernel
}  // namespace noncentral
}  // namespace openMVG


namespace openMVG{
namespace SfMRobust{

  using namespace openMVG::matching;
  using namespace opengv;
  using namespace std;

/**
 * @brief Estimate the relative rig pose from normalized point matches.
 *
 * @param[in] b1 bearing vectors of rig one
 * @param[in] b2 bearing vectors of rig two
 * @param[in] scIdOne subcamera id of each bearing vector of rig one
 * @param[in] scIdTwo subcamera id of each bearing vector of rig two
 * @param[in] rigOffsets center of cameras in rig rig referential
 * @param[in] rigRotations rotation matrices of subcameras
 * @param[out] transformation_t relative pose of second rig (R^T and C)
 * @param[out] pvec_inliers inliers indices (can be empty)
 * @param[out] errorMax upper bound of the reprojection error of the found solution
 * @param[in] precision upper bound of the desired solution
 */
bool robustRigPose(
  const bearingVectors_t & b1,
  const bearingVectors_t & b2,
  const std::vector<int> & scIdOne,
  const std::vector<int> & scIdTwo,
  const translations_t & rigOffsets,
  const rotations_t & rigRotations,
  transformation_t * relativePose,
  std::vector<size_t> * pvec_inliers,
  double * errorMax,
  double precision);

} // namespace SfMRobust
} // namespace openMVG
