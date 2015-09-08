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

#ifndef OPENMVGRIG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
#define OPENMVGRIG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/graph/graph.hpp"

using namespace openMVG::rotation_averaging;

namespace openMVG{
namespace sfm{

class GlobalSfMRig_Translation_AveragingSolver
{
  RelativeInfo_Vec _vec_initialRijTijEstimates;

public:

  /// Use features in normalized camera frames
  bool Run(
    ETranslationAveragingMethod eTranslationAveragingMethod,
    SfM_Data & sfm_data,
    const Features_Provider * normalized_features_provider,
    const Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const openMVG::rotation_averaging::RelativeRotations & vec_relatives_R,
    matching::PairWiseMatches & tripletWise_matches
  );

private:
  bool Translation_averaging(
    ETranslationAveragingMethod eTranslationAveragingMethod,
    SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR);

  void Compute_translations(
    const SfM_Data & sfm_data,
    const Features_Provider * normalized_features_provider,
    const Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const openMVG::rotation_averaging::RelativeRotations & vec_relatives_R,
    matching::PairWiseMatches &tripletWise_matches);

  //-- Compute the relative translations on the rotations graph.
  // Compute relative translations by using triplets of poses.
  // Use an edge coverage algorithm to reduce the graph covering complexity
  // Complexity: sub-linear in term of edges count.
  void ComputePutativeTranslation_EdgesCoverage(
    const SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const openMVG::rotation_averaging::RelativeRotations & vec_relatives_R,
    const Features_Provider * normalized_features_provider,
    const Matches_Provider * matches_provider,
    RelativeInfo_Vec & vec_initialEstimates,
    matching::PairWiseMatches & newpairMatches);

  // Robust estimation and refinement of a translation and 3D points of an image triplets.
  bool Estimate_T_triplet(
    const SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const Features_Provider * normalized_features_provider,
    const Matches_Provider * matches_provider,
    const graph::Triplet & poses_id,
    std::vector<Vec3> & vec_tis,
    double & dPrecision, // UpperBound of the precision found by the AContrario estimator
    std::vector<size_t> & vec_inliers,
    openMVG::tracks::STLMAPTracks & rig_tracks,
    const double ThresholdUpperBound, //Threshold used for the trifocal tensor estimation solver used in AContrario Ransac
    const std::string & sOutDirectory) const;

  // triplet rotation rejection based on angular threshold
  void TripletRotationRejection(
    const double max_angular_error,
    std::vector< graph::Triplet > & vec_triplets,
    RelativeRotations & relativeRotations) const;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
