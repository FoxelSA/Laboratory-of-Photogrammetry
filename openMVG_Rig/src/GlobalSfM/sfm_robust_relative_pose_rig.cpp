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

#include  "./sfm_robust_relative_pose_rig.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include "openMVG/numeric/poly.h"

#include <set>


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

/*************************
 *
 *  Definition of solver used for non central relative pose estimation
 *
 ************************
*/
namespace openMVG {
namespace noncentral {
namespace kernel {

using namespace std;
using namespace opengv;

/**
 * Six point solver for non central camera system,
 * // [1] "Solutions to minimal generalized relative pose problems".
 * // authors: Stewenius, H., Nister, D., Oskarsson, M., & Astrom K,
 * // Date: 2005:
 * // Conference: Workshop on omnidirectional vision 2005.
 */
void SixPointSolver::Solve(
                  relative_pose::NoncentralRelativeAdapter & adapter,
                  std::vector<transformation_t> * models,
                  const std::vector<size_t> &indices)
{

  // convert size_t to int for opengv call
  std::vector<int> idx(indices.begin(), indices.end());

  // create non central relative sac problem
  sac_problems::relative_pose::NoncentralRelativePoseSacProblem
            problem(adapter,
                    sac_problems::relative_pose::NoncentralRelativePoseSacProblem::SIXPT,
                    false);

  // solve pose problem
  transformation_t relativePose;
  problem.computeModelCoefficients(idx, relativePose);

  models->push_back(relativePose);
}

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

void GePointSolver::Solve(
          relative_pose::NoncentralRelativeAdapter & adapter,
          std::vector<transformation_t> * models,
          const std::vector<size_t> &indices)
{

  // convert size_t to int for opengv call
  std::vector<int> idx(indices.begin(), indices.end());

  // create non central relative sac problem
  sac_problems::relative_pose::NoncentralRelativePoseSacProblem
            problem(adapter,
                    sac_problems::relative_pose::NoncentralRelativePoseSacProblem::GE,
                    false);

  // solve pose problem
  transformation_t relativePose;
  problem.computeModelCoefficients(idx, relativePose);

  models->push_back(relativePose);
}

}  // namespace kernel
}  // namespace noncentral
}  // namespace openMVG


/*******************************************************************************
 *
 *  Robust rig pose
 *
 *******************************************************************************
 */

namespace openMVG{
namespace SfMRobust{

  using namespace openMVG::matching;
  using namespace openMVG::robust;
  using namespace opengv;

static const size_t ACRANSAC_ITER = 4096;


/**
 *  Robust rig pose estimation
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
  double precision = std::numeric_limits<double>::infinity())
{
  assert(pvec_inliers != NULL);

  // Use the Generalized eigenvalue solver to solve pose problem
  typedef openMVG::noncentral::kernel::GePointSolver SolverType;

  // Define the AContrario adaptor
  typedef ACKernelAdaptorRigPose<
      SolverType,
      openMVG::noncentral::kernel::RigAngularError,
      transformation_t>
      KernelType;

  KernelType kernel(b1, b2, scIdOne, scIdTwo, rigOffsets, rigRotations);

  // Robustly estimation of the Essential matrix and it's precision
  std::pair<double,double> acRansacOut = ACRANSAC(kernel, *pvec_inliers,
    ACRANSAC_ITER, relativePose, precision, false );
  *errorMax = acRansacOut.first;

  return ( pvec_inliers->size() > 2.5 * SolverType::MINIMUM_SAMPLES * rigOffsets.size() );
}

} // namespace SfMRobust
} // namespace openMVG
