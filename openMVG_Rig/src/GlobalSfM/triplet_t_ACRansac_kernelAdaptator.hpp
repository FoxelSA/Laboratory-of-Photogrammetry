/*
 * triplet_t_ACRansac_kernelAdaptator
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
#define USE_L_INFINITY_TRANSLATION 0

#include "./triplet_t_solver_rig.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/conditioning.hpp"

#include "openMVG/linearProgramming/linearProgramming.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"

namespace openMVG {
  namespace trifocal {
    namespace kernel {

      /// A trifocal tensor seen as 3 projective cameras
      struct rigTrackTrifocalTensorModel {
        Mat3  R1, R2, R3;
        Vec3  t1, t2, t3;

        static double Error(const rigTrackTrifocalTensorModel & t,
        std::vector < std::vector <double> > pointInfo,
        const std::vector<Mat3> & rigRotation, // rig subcamera rotation
        const std::vector<Vec3> & rigOffsets ) // rig subcamera translation
      {
        // Triangulate and return the reprojection error
        Triangulation triangulationObj;
        std::vector < std::pair <Mat34, Vec2> >  views;

        for( size_t i = 0 ; i < pointInfo.size() ; ++ i)
        {
          // extract sucamera rotations and translation
          const size_t I = (size_t) pointInfo[i][2]; // intrinsic id
          const Mat3  RI = rigRotation[I];  const Vec3 tI = -RI * rigOffsets[I];

          // compute projection matrix
          Mat34 P ;
          switch( (int) pointInfo[i][3] )  // pose Id
          {
            // if first pose
            case 0 :
              P = HStack(RI * t.R1, RI * t.t1 + tI);
              break;

            // if second pose
            case 1 :
              P = HStack(RI * t.R2, RI * t.t2 + tI);
              break;

            // if third pose
            case 2 :
              P = HStack(RI * t.R3, RI * t.t3 + tI);
              break;
          };

          //compute projection matrices
          Vec2 pt;
          pt << pointInfo[i][0], pointInfo[i][1];

          views.push_back(std::make_pair ( P, pt ) );
        }

        // update triangulation object
        for( size_t i = 0 ; i < views.size(); ++i )
          triangulationObj.add ( views[i].first, views[i].second );

        const Vec3 X = triangulationObj.compute();

        //- Return error
        double max_error = 0.0;

        // update triangulation object
        for( size_t i = 0 ; i < views.size(); ++i )
        {
          const Mat34 P = views[i].first;
          const Vec2 pt = views[i].second;
          max_error = std::max( (Project(P, X) - pt ).norm(), max_error );
        }

        //- Return max error as a test
        return max_error;
      }
    };

  }  // namespace kernel
}  // namespace trifocal
}  // namespace openMVG


namespace openMVG{

  using namespace openMVG::trifocal::kernel;

  struct rigTrackTisXisTrifocalSolver {
    enum { MINIMUM_SAMPLES = 4 };
    enum { MAX_MODELS = 1 };
    // Solve the computation of the tensor.
    static void Solve(
    const std::vector< std::vector < std::vector <double> > > pt,
    const std::vector<Mat3> & rigRotation, // rotation of subcameras
    const std::vector<Vec3> & rigOffsets, // optical center of rig subcameras in rig referential frame
    const std::vector<Mat3> & vec_KR, std::vector<rigTrackTrifocalTensorModel> *P,
    const double ThresholdUpperBound)
  {
    //Build the megaMatMatrix
    int n_obs = 0;
    for(int i=0; i < pt.size() ; ++i )
      n_obs += pt[i].size();

    Mat megaMat(5, n_obs);
    {
      size_t cpt = 0;
      for (size_t i = 0; i  < pt.size(); ++i)
      {
        for(size_t j = 0; j < pt[i].size(); ++j)
        {
          megaMat.col(cpt) << pt[i][j][0], pt[i][j][1], i, pt[i][j][2], pt[i][j][3]; // feature x and y, 3d point index, rig subcam index and rig index
          ++cpt;
        }
      }
    }

#if USE_L_INFINITY_TRANSLATION
      //-- Solve the LInfinity translation and structure from Rotation and points data.
      std::vector<double> vec_solution((3 + MINIMUM_SAMPLES)*3);

      using namespace openMVG::lInfinityCV;

      #ifdef OPENMVG_HAVE_MOSEK
      MOSEK_SolveWrapper LPsolver(static_cast<int>(vec_solution.size()));
      #else
      OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));
      #endif

      Rig_Translation_Structure_L1_ConstraintBuilder cstBuilder(vec_KR, megaMat, rigRotation, rigOffsets);
      double gamma;
      if (BisectionLP<Rig_Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
        LPsolver,
        cstBuilder,
        &vec_solution,
        ThresholdUpperBound,//admissibleResidual,
        0.0, 1e-8, 2, &gamma, false))
      {
        std::vector<Vec3> vec_tis(3);
        vec_tis[0] = Vec3(vec_solution[0], vec_solution[1], vec_solution[2]);
        vec_tis[1] = Vec3(vec_solution[3], vec_solution[4], vec_solution[5]);
        vec_tis[2] = Vec3(vec_solution[6], vec_solution[7], vec_solution[8]);

        rigTrackTrifocalTensorModel PTemp;
        PTemp.R1 = vec_KR[0]; PTemp.t1 = vec_tis[0];
        PTemp.R2 = vec_KR[1]; PTemp.t2 = vec_tis[1];
        PTemp.R3 = vec_KR[2]; PTemp.t3 = vec_tis[2];

        P->push_back(PTemp);
      }
#else
      //-- Solve the LInfinity translation and structure from Rotation and points data.
      std::vector<double> vec_solution((9 + n_obs));

      using namespace openMVG::lInfinityCV;

      #ifdef OPENMVG_HAVE_MOSEK
      MOSEK_SolveWrapper LPsolver(static_cast<int>(vec_solution.size()));
      #else
      OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));
      #endif

      Rig_Center_Structure_L1_ConstraintBuilder cstBuilder(vec_KR, megaMat, rigRotation, rigOffsets);
      double gamma;
      if (BisectionLP<Rig_Center_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
        LPsolver,
        cstBuilder,
        &vec_solution,
        ThresholdUpperBound,//admissibleResidual,
        0.0, 1e-8, 2, &gamma, false))
      {
        std::vector<Vec3> vec_tis(3);
        vec_tis[0] = Vec3(vec_solution[0], vec_solution[1], vec_solution[2]);
        vec_tis[1] = Vec3(vec_solution[3], vec_solution[4], vec_solution[5]);
        vec_tis[2] = Vec3(vec_solution[6], vec_solution[7], vec_solution[8]);

        rigTrackTrifocalTensorModel PTemp;
        PTemp.R1 = vec_KR[0]; PTemp.t1 = vec_tis[0];
        PTemp.R2 = vec_KR[1]; PTemp.t2 = vec_tis[1];
        PTemp.R3 = vec_KR[2]; PTemp.t3 = vec_tis[2];

        P->push_back(PTemp);
      }
#endif
    }

    // Compute the residual of reprojections
    static double Error(const rigTrackTrifocalTensorModel & Tensor, const std::vector < std::vector <double> > featInfo,
    const std::vector<Mat3> & rigRotation, const std::vector<Vec3> & rigOffsets)
  {
    return rigTrackTrifocalTensorModel::Error(Tensor, featInfo, rigRotation, rigOffsets);
  }
};

template <typename SolverArg,
typename ErrorArg,
typename ModelArg>
class rig_TrackTrifocalKernel_ACRansac_N_tisXis
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;


  rig_TrackTrifocalKernel_ACRansac_N_tisXis(const std::vector< std::vector < std::vector <double> > > pt,
  const std::vector<Mat3> & vec_KRi,
  const std::vector<Mat3> & rigRotation,
  const std::vector<Vec3> & rigOffsets,
  const double ThresholdUpperBound,
  const std::pair<IndexT, IndexT> image_dimension)
  : pt_(pt), vec_KR_(vec_KRi),
  vec_rigRotation_(rigRotation),
  vec_rigOffset_(rigOffsets),
  ThresholdUpperBound_(ThresholdUpperBound),
  logalpha0_(log10(M_PI/(image_dimension.first*image_dimension.second)))
{
  //initialize normalized coordinates
  // Normalize points by inverse(K)
}

enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
enum { MAX_MODELS = Solver::MAX_MODELS };

void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {

  std::vector < std::vector < std::vector <double > > > pt_sampled;

  for( int i=0; i < samples.size(); ++i )
    pt_sampled.push_back( pt_[samples[i]] );

    // Create a model from the points
    Solver::Solve( pt_sampled,
    vec_rigRotation_, vec_rigOffset_,
    vec_KR_, models, ThresholdUpperBound_);
  }

  void Errors(const Model &model, std::vector<double> & vec_errors) const {
    for (size_t sample = 0; sample < pt_.size(); ++sample)
    {
      vec_errors[sample] = ErrorArg::Error(model, pt_[sample], vec_rigRotation_, vec_rigOffset_);
    }
  }

  double Error(size_t sample, const Model &model) const {
    return ErrorArg::Error(model, pt_[sample], vec_rigRotation_, vec_rigOffset_);
  }

  size_t NumSamples() const {
    return pt_.size();
  }

  void Unnormalize(Model * model) const {
    // Unnormalize model is not necessary since K is considered as Identity.
  }

  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const { return val;}

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 1.0;}

private:
  const std::vector < std::vector < std::vector <double > > > pt_;
  const double logalpha0_;
  const double ThresholdUpperBound_;
  const std::vector<Mat3> vec_KR_;
  const std::vector<Mat3> vec_rigRotation_;
  const std::vector<Vec3> vec_rigOffset_;

};
} // namespace openMVG
