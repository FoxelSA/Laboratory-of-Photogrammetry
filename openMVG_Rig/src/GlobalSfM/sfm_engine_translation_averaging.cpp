/*
 * sfm_engine_translation_averaging
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

#include "./sfm_engine_translation_averaging.hpp"

#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/sfm/pipelines/global/mutexSet.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/linearProgramming/linearProgramming.hpp"

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/conditioning.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/sfm/pipelines/global/triplet_t_ACRansac_kernelAdaptator.hpp"


#include "third_party/histogram/histogram.hpp"
#include "third_party/progress/progress.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;

/// Use features in normalized camera frames
bool GlobalSfMRig_Translation_AveragingSolver::Run(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  SfM_Data & sfm_data,
  const Features_Provider * normalized_features_provider,
  const Matches_Provider * matches_provider,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Compute the relative translations and save them to vec_initialRijTijEstimates:
  Compute_translations(
    sfm_data,
    normalized_features_provider,
    matches_provider,
    map_globalR,
    tripletWise_matches);

  // TODO: Keep the largest Biedge connected component graph of relative translations

  // Compute the global translations
  return Translation_averaging(
    eTranslationAveragingMethod,
    sfm_data,
    map_globalR);
}

bool GlobalSfMRig_Translation_AveragingSolver::Translation_averaging(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR)
{
  return true;
}

void GlobalSfMRig_Translation_AveragingSolver::Compute_translations(
  const SfM_Data & sfm_data,
  const Features_Provider * normalized_features_provider,
  const Matches_Provider * matches_provider,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  matching::PairWiseMatches &tripletWise_matches)
{
  std::cout << "\n-------------------------------" << "\n"
    << " Relative translations computation: " << "\n"
    << "-------------------------------" << std::endl;

  // Compute relative translations over the graph of global rotations
  //  thanks to an edge coverage algorithm


  ComputePutativeTranslation_EdgesCoverage(
    sfm_data,
    map_globalR,
    normalized_features_provider,
    matches_provider,
    vec_initialRijTijEstimates,
    tripletWise_matches);
}

//-- Perform a trifocal estimation of the graph contain in vec_triplets with an
// edge coverage algorithm. It's complexity is sub-linear in term of edges count.
void GlobalSfMRig_Translation_AveragingSolver::ComputePutativeTranslation_EdgesCoverage(
  const SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  const Features_Provider * normalized_features_provider,
  const Matches_Provider * matches_provider,
  RelativeInfo_Vec & vec_initialEstimates,
  matching::PairWiseMatches & newpairMatches)
{
  openMVG::system::Timer timerLP_triplet;

  //--
  // Compute the relative translations using triplets of rotations over the rotation graph.
  //--
  //
  // 1. List plausible triplets over the global rotation pose graph Ids.
  //   - list all edges that have support in the rotation pose graph
  //
  Pair_Set matches_pair_belong_to_pose_ids;
  std::set<IndexT> set_pose_ids;
  std::transform(map_globalR.begin(), map_globalR.end(),
    std::inserter(set_pose_ids, set_pose_ids.begin()), stl::RetrieveKey());
  // List shared correspondences (pairs) between poses
  for (const auto & match_iterator : matches_provider->_pairWise_matches)
  {
    const Pair pair = match_iterator.first;
    const View * v1 = sfm_data.GetViews().at(pair.first).get();
    const View * v2 = sfm_data.GetViews().at(pair.second).get();
    // Consider iff the pair is supported by the rotation graph
    if (v1->id_pose != v2->id_pose &&
        set_pose_ids.count(v1->id_pose) && set_pose_ids.count(v2->id_pose))
    {
      matches_pair_belong_to_pose_ids.insert(
        std::move(std::make_pair(v1->id_pose, v2->id_pose)));
    }
  }
  // List putative triplets (from global rotations Ids)
  const std::vector< graph::Triplet > vec_triplets =
    graph::tripletListing(matches_pair_belong_to_pose_ids);
  std::cout << "#Triplets: " << vec_triplets.size() << std::endl;


  {
    // Compute triplets of translations
    // Avoid to cover each edge of the graph by using an edge coverage algorithm
    // An estimated triplets of translation mark three edges as estimated.
  }

  const double timeLP_triplet = timerLP_triplet.elapsed();
  std::cout << "TRIPLET COVERAGE TIMING: " << timeLP_triplet << " seconds" << std::endl;

  std::cout << "-------------------------------" << "\n"
      << "-- #Effective translations estimates: " << vec_initialRijTijEstimates.size()/3
      << " from " << vec_triplets.size() << " triplets.\n"
      << "-- resulting in " <<vec_initialRijTijEstimates.size() << " translation estimation.\n"
      << "-- timing to obtain the relative translations: " << timeLP_triplet << " seconds.\n"
      << "-------------------------------" << std::endl;
}

// Robust estimation and refinement of a translation and 3D points of an image triplets.
bool GlobalSfMRig_Translation_AveragingSolver::Estimate_T_triplet(
  const SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  const Features_Provider * normalized_features_provider,
  const openMVG::tracks::STLMAPTracks & map_tracksCommon,
  const graph::Triplet & poses_id,
  std::vector<Vec3> & vec_tis,
  double & dPrecision, // UpperBound of the precision found by the AContrario estimator
  std::vector<size_t> & vec_inliers,
  const double ThresholdUpperBound, //Threshold used for the trifocal tensor estimation solver used in AContrario Ransac
  const std::string & sOutDirectory) const
{
  return false;
}

} // namespace sfm
} // namespace openMVG
