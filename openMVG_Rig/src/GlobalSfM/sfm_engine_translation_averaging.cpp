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
#include "./triplet_t_ACRansac_kernelAdaptator.hpp"

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
  Pair_Set rotation_pose_id_graph;
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
      rotation_pose_id_graph.insert(
        std::move(std::make_pair(v1->id_pose, v2->id_pose)));
    }
  }
  // List putative triplets (from global rotations Ids)
  const std::vector< graph::Triplet > vec_triplets =
    graph::tripletListing(rotation_pose_id_graph);
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
  const Matches_Provider * matches_provider,
  const graph::Triplet & poses_id,
  std::vector<Vec3> & vec_tis,
  double & dPrecision, // UpperBound of the precision found by the AContrario estimator
  std::vector<size_t> & vec_inliers,
  const double ThresholdUpperBound, //Threshold used for the trifocal tensor estimation solver used in AContrario Ransac
  const std::string & sOutDirectory) const
{
  // namespaces
  using namespace linearProgramming;
  using namespace lInfinityCV;
  using namespace openMVG::trifocal;
  using namespace openMVG::trifocal::kernel;
  using namespace openMVG::tracks;

  // poses index initialization
  const size_t I = poses_id.i, J = poses_id.j , K = poses_id.k;

  // List matches that belong to the triplet of poses
  PairWiseMatches map_triplet_matches;
  std::set<IndexT> set_pose_ids;
  set_pose_ids.insert(poses_id.i);
  set_pose_ids.insert(poses_id.j);
  set_pose_ids.insert(poses_id.k);
  // List shared correspondences (pairs) between poses
  for (const auto & match_iterator : matches_provider->_pairWise_matches)
  {
    const Pair pair = match_iterator.first;
    const View * v1 = sfm_data.GetViews().at(pair.first).get();
    const View * v2 = sfm_data.GetViews().at(pair.second).get();
    // Consider iff the pair is supported by the triplet & rotation graph
    const bool b_different_pose_id = v1->id_pose != v2->id_pose;
    const int covered_pose =
      set_pose_ids.count(v1->id_pose) +
      set_pose_ids.count(v2->id_pose);
    // Different pose Id and the current edge cover the triplet edge
    if (b_different_pose_id && covered_pose == 2 )
    {
      map_triplet_matches.insert(
        std::make_pair( pair, std::move(match_iterator.second)) );
    }
  }

  openMVG::tracks::STLMAPTracks rig_tracks;

  openMVG::tracks::TracksBuilder tracksBuilder;
  tracksBuilder.Build(map_triplet_matches);
  tracksBuilder.Filter(3);
  tracksBuilder.ExportToSTL(rig_tracks);

  // Evaluate the triplet for relative translation computation:
  // 1. From tracks initialize the structure observation
  // 2. Setup the known parameters (camera intrinsics K + subposes) + global rotations
  // 3. Solve the unknown: relative translations

  // Get rotations:
  std::vector<Mat3> vec_global_KR_Triplet;
  vec_global_KR_Triplet.push_back(map_globalR.at(I));
  vec_global_KR_Triplet.push_back(map_globalR.at(J));
  vec_global_KR_Triplet.push_back(map_globalR.at(K));

  // check that there is enough correspondances to evaluate model
  const size_t  rigSize = sfm_data.GetIntrinsics().size();
  if( rig_tracks.size() < 50 * rigSize )
    return false ;

  // initialize rig structure for relative translation estimation
  std::vector<Vec3>  rigOffsets;
  std::vector<Mat3>  rigRotations;
  double             minFocal=1.0e10;

  // Update rig structure from OpenMVG data to OpenGV convention
  for (const auto & intrinsicVal : sfm_data.GetIntrinsics())
  {
    const cameras::IntrinsicBase * intrinsicPtr = intrinsicVal.second.get();
    if ( intrinsicPtr->getType() == cameras::PINHOLE_RIG_CAMERA )
    {
      // retrieve information from shared pointer
      const cameras::Rig_Pinhole_Intrinsic * rig_intrinsicPtr = dynamic_cast< const cameras::Rig_Pinhole_Intrinsic * > (intrinsicPtr);
      const geometry::Pose3 sub_pose = rig_intrinsicPtr->get_subpose();
      const double focal = rig_intrinsicPtr->focal();

      // update rig stucture
      const IndexT index = intrinsicVal.first;
      rigOffsets[index]   = sub_pose.center();
      rigRotations[index] = sub_pose.rotation();

      minFocal = std::min( minFocal , focal );
    }
  }

  // initialize rigId map
  std::map  < size_t, size_t > map_rigIdToTripletId;
  map_rigIdToTripletId[I] = 0; map_rigIdToTripletId[J] = 1; map_rigIdToTripletId[K] = 2;

  // initialize data for model evaluation
  std::vector < std::vector < std::vector < double > > > featsAndRigIdPerTrack;
  std::map  <size_t, size_t>  sampleToTrackId;
  size_t cpt = 0;

  // List the tracks to associate a pair of bearing vector to a track Id
  std::map < size_t, size_t >  map_bearingIdToTrackId;
  for (STLMAPTracks::const_iterator iterTracks = rig_tracks.begin();
    iterTracks != rig_tracks.end(); ++iterTracks, ++cpt)
  {
    const submapTrack & track = iterTracks->second;
    for (submapTrack::const_iterator iterTrack_I = track.begin();
      iterTrack_I != track.end(); ++iterTrack_I)
    {
      const size_t I  = iterTrack_I->first;
      const size_t feat_I = iterTrack_I->second;
      submapTrack::const_iterator iterTrack_J = iterTrack_I;
      std::advance(iterTrack_J, 1);

      for (iterTrack_J; iterTrack_J != track.end(); ++iterTrack_J)
      {
        const size_t J  = iterTrack_J->first;
        const size_t feat_J = iterTrack_J->second;
        submapTrack::const_iterator iterTrack_K = iterTrack_J;
        std::advance(iterTrack_K, 1);

        for (iterTrack_K; iterTrack_K != track.end(); ++iterTrack_K)
        {
          const size_t K  = iterTrack_K->first;
          const size_t feat_K = iterTrack_K->second;

          // initialize view structure
          const View * view_I = sfm_data.views.at(I).get();
          const View * view_J = sfm_data.views.at(J).get();
          const View * view_K = sfm_data.views.at(K).get();

          // initialize intrinsic group of cameras I and J
          const IndexT intrinsic_index_I = view_I->id_intrinsic;
          const IndexT intrinsic_index_J = view_J->id_intrinsic;
          const IndexT intrinsic_index_K = view_K->id_intrinsic;

          // extract normalized keypoints coordinates
          Vec3   bearing_I, bearing_J, bearing_K;
          bearing_I << normalized_features_provider->feats_per_view.at(I).at(feat_I).coords().cast<double>(), 1.0;
          bearing_J << normalized_features_provider->feats_per_view.at(J).at(feat_J).coords().cast<double>(), 1.0;
          bearing_K << normalized_features_provider->feats_per_view.at(K).at(feat_K).coords().cast<double>(), 1.0;

          // initialize relative translation data container
          std::vector<double>   feat_cam_I(4);
          std::vector<double>   feat_cam_J(4);
          std::vector<double>   feat_cam_K(4);

          // fill data container
          feat_cam_I[0] = bearing_I[0];  feat_cam_I[1] = bearing_I[1]; feat_cam_I[2] = intrinsic_index_I; feat_cam_I[3] = map_rigIdToTripletId.at(I);
          feat_cam_J[0] = bearing_J[0];  feat_cam_J[1] = bearing_J[1]; feat_cam_J[2] = intrinsic_index_J; feat_cam_J[3] = map_rigIdToTripletId.at(J);
          feat_cam_K[0] = bearing_K[0];  feat_cam_K[1] = bearing_K[1]; feat_cam_K[2] = intrinsic_index_K; feat_cam_K[3] = map_rigIdToTripletId.at(K);

          // export it
          std::vector < std::vector < double > >   tmp;
          tmp.push_back(feat_cam_I);
          tmp.push_back(feat_cam_J);
          tmp.push_back(feat_cam_K);
          featsAndRigIdPerTrack.push_back( tmp );
          sampleToTrackId[ sampleToTrackId.size() ] = cpt;
        }
      }
    }
  }
  // set thresholds for relative translation estimation
  const size_t  ORSA_ITER = 1024;             // max number of iterations of AC-RANSAC

  // compute model
  typedef  rigTrackTisXisTrifocalSolver  SolverType;

  typedef rig_TrackTrifocalKernel_ACRansac_N_tisXis<
    rigTrackTisXisTrifocalSolver,
    rigTrackTisXisTrifocalSolver,
    rigTrackTrifocalTensorModel> KernelType;
  KernelType kernel(featsAndRigIdPerTrack, vec_global_KR_Triplet, rigRotations, rigOffsets, ThresholdUpperBound);

  rigTrackTrifocalTensorModel T;
  std::pair<double,double> acStat = robust::ACRANSAC(kernel, vec_inliers, ORSA_ITER, &T, dPrecision, false );
  dPrecision = acStat.first;

  //-- Export data in order to have an idea of the precision of the estimates
  vec_tis.resize(3);
  vec_tis[0] = T.t1;
  vec_tis[1] = T.t2;
  vec_tis[2] = T.t3;

  // update inlier list
  std::set <size_t>  inliers_tracks;
  for( size_t i = 0 ; i < vec_inliers.size() ; ++i )
      inliers_tracks.insert( sampleToTrackId[vec_inliers[i]] );

  std::vector <size_t>  inliers;
  for( size_t i = 0 ; i < rig_tracks.size() ; ++i )
     if( inliers_tracks.find(i) !=  inliers_tracks.end() )
          inliers.push_back( i );

  vec_inliers.swap( inliers );


  // if there is more than 2/3 of inliers, keep mondel
  const bool bTest =  ( vec_inliers.size() > 0.66 * rig_tracks.size() ) ;

  if (!bTest)
  {
    std::cout << "Triplet rejected : AC: " << dPrecision
      << " inliers count " << inliers_tracks.size()
      << " total putative " << featsAndRigIdPerTrack.size() << std::endl;
  }

  return bTest ;

}

} // namespace sfm
} // namespace openMVG
