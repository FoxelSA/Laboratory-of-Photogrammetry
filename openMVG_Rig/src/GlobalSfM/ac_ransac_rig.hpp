// Copyright (c) 2012, 2013 Lionel MOISAN.
// Copyright (c) 2012, 2013 Pascal MONASSE.
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
//-------------------
// Generic implementation of ACRANSAC
//-------------------
// The A contrario parametrization have been first explained in [1] and
//  later extended to generic model estimation in [2] (with a demonstration for
//  the homography) and extended and use at large scale for Structure from
//  Motion in [3].
//
//--
//  [1] Lionel Moisan, Berenger Stival,
//  A probalistic criterion to detect rigid point matches between
//  two images and estimate the fundamental matrix.
//  IJCV 04.
//--
//  [2] Lionel Moisan, Pierre Moulon, Pascal Monasse.
//  Automatic Homographic Registration of a Pair of Images,
//    with A Contrario Elimination of Outliers
//  Image Processing On Line (IPOL), 2012.
//  http://dx.doi.org/10.5201/ipol.2012.mmm-oh
//--
//  [3] Pierre Moulon, Pascal Monasse and Renaud Marlet.
//  Adaptive Structure from Motion with a contrario mode estimation.
//  In 11th Asian Conference on Computer Vision (ACCV 2012)
//--

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/robust_estimation/rand_sampling.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
#include <limits>
#include <iostream>

namespace openMVG {
namespace robust{
/**
 * @brief ACRANSAC routine (ErrorThreshold, NFA)
 *
 * @param[in] kernel model and metric object
 * @param[out] vec_inliers points that fit the estimated model
 * @param[in] nIter maximum number of consecutive iterations
 * @param[out] model returned model if found
 * @param[in] precision upper bound of the precision (squared error)
 * @param[in] bVerbose display console log
 *
 * @return (errorMax, minNFA)
 */
template<typename Kernel>
std::pair<double, double> ACRANSAC_RIG(const Kernel &kernel,
  std::vector<size_t> & vec_inliers,
  size_t nIter = 1024,
  typename Kernel::Model * model = NULL,
  double precision = std::numeric_limits<double>::infinity(),
  bool bVerbose = false )
{
  const size_t sizeSample = Kernel::MINIMUM_SAMPLES;
  const size_t nData = kernel.NumSamples();
  if(nData <= (size_t)sizeSample)
    return std::make_pair(0.0,0.0);

  // Output parameters
  double minNFA = std::numeric_limits<double>::infinity();
  double errorMax = std::numeric_limits<double>::infinity();

  vec_inliers.clear();

  const double maxThreshold = (precision==std::numeric_limits<double>::infinity()) ?
    std::numeric_limits<double>::infinity() :
    precision * kernel.normalizer2()(0,0) * kernel.normalizer2()(0,0);

  std::vector<ErrorIndex> vec_residuals(nData); // [residual,index]
  std::vector<double> vec_residuals_(nData);
  std::vector<size_t> vec_sample(sizeSample); // Sample indices

  // Possible sampling indices (could change in the optimization phase)
  std::vector<size_t> vec_index(nData);
  for (size_t i = 0; i < nData; ++i)
    vec_index[i] = i;

  // Precompute log combi
  double loge0 = log10((double)Kernel::MAX_MODELS * (nData-sizeSample));
  std::vector<float> vec_logc_n, vec_logc_k;
  makelogcombi_n(nData, vec_logc_n);
  makelogcombi_k(sizeSample, nData, vec_logc_k);

  // Reserve 10% of iterations for focused sampling
  const size_t   given_N = nIter ;
  size_t nIterReserve = nIter/10;
  nIter -= nIterReserve;

  // Main estimation loop.
  for (size_t iter=0; iter < nIter; ++iter) {
    UniformSample(sizeSample, vec_index, &vec_sample); // Get random sample

    std::vector<typename Kernel::Model> vec_models; // Up to max_models solutions
    kernel.Fit(vec_sample, &vec_models);

    // Evaluate models
    bool better = false;
    for (size_t k = 0; k < vec_models.size(); ++k)  {
      // Residuals computation and ordering
      kernel.Errors(vec_models[k], vec_residuals_);
      for (size_t i = 0; i < nData; ++i)  {
        const double error = vec_residuals_[i];
        vec_residuals[i] = ErrorIndex(error, i);
      }
      std::sort(vec_residuals.begin(), vec_residuals.end());

      // Most meaningful discrimination inliers/outliers
      const ErrorIndex best = bestNFA(
        sizeSample,
        kernel.logalpha0(),
        vec_residuals,
        loge0,
        maxThreshold,
        vec_logc_n,
        vec_logc_k,
        kernel.multError());

      if (best.first < minNFA  && best.second > vec_inliers.size() )  {
        // A better model was found
        better = true;
        minNFA = best.first;
        vec_inliers.resize(best.second);
        for (size_t i=0; i<best.second; ++i)
          vec_inliers[i] = vec_residuals[i].second;
        errorMax = vec_residuals[best.second-1].first; // Error threshold
        if(model) *model = vec_models[k];

        if(bVerbose)  {
          std::cout << "  nfa=" << minNFA
            << " inliers=" << best.second
            << " precisionNormalized=" << errorMax
            << " precision=" << kernel.unormalizeError(errorMax)
            << " (iter=" << iter;
          std::cout << ",sample=";
          std::copy(vec_sample.begin(), vec_sample.end(),
            std::ostream_iterator<size_t>(std::cout, ","));
          std::cout << ")" <<std::endl;
        }
      }
    }

    // ACRANSAC optimization: draw samples among best set of inliers so far
    if((better && minNFA<0) || (iter+1==nIter && nIterReserve)) {
      if(vec_inliers.empty()) { // No model found at all so far
        nIter++; // Continue to look for any model, even not meaningful
        nIterReserve--;
      } else
      {
        if( vec_inliers.size() > 0.85 * nData )
        {
          // ACRANSAC optimization: draw samples among best set of inliers so far
          vec_index = vec_inliers;
          if(nIterReserve) {
              nIter = iter+1+nIterReserve;
              nIterReserve=0;
          }
        }
      }
    }
  }

  if(minNFA >= 0)
    vec_inliers.clear();

  if (!vec_inliers.empty())
  {
    if (model)
      kernel.Unnormalize(model);
    errorMax = kernel.unormalizeError(errorMax);
  }

  return std::make_pair(errorMax, minNFA);
}

} // namespace robust
} // namespace openMVG
