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
        std::make_pair(v1->id_pose, v2->id_pose));
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

    //-- precompute the number of track per triplet:
    Hash_Map<IndexT, IndexT> map_tracksPerTriplets;
    #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < (int)vec_triplets.size(); ++i)
    {
      // List matches that belong to the triplet of poses
      const graph::Triplet & triplet = vec_triplets[i];
      PairWiseMatches map_triplet_matches;
      std::set<IndexT> set_pose_ids;
      set_pose_ids.insert(triplet.i);
      set_pose_ids.insert(triplet.j);
      set_pose_ids.insert(triplet.k);
      // List shared correspondences (pairs) between the triplet poses
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
          map_triplet_matches.insert( match_iterator );
        }
      }
      // Compute tracks:
      {
        openMVG::tracks::TracksBuilder tracksBuilder;
        tracksBuilder.Build(map_triplet_matches);
        tracksBuilder.Filter(3);
        #ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
        #endif
        map_tracksPerTriplets[i] = tracksBuilder.NbTracks(); //count the # of matches in the UF tree
      }
    }

    typedef Pair myEdge;

    //-- List all edges
    std::set<myEdge > set_edges;

    for (size_t i = 0; i < vec_triplets.size(); ++i)
    {
      const graph::Triplet & triplet = vec_triplets[i];
      const IndexT I = triplet.i, J = triplet.j , K = triplet.k;
      // Add three edges
      set_edges.insert(std::make_pair(I,J));
      set_edges.insert(std::make_pair(I,K));
      set_edges.insert(std::make_pair(J,K));
    }
    // Move set to a vector
    std::vector<myEdge > vec_edges(std::begin(set_edges), std::end(set_edges));
    std::set<myEdge >().swap(set_edges); // release memory

    openMVG::sfm::MutexSet<myEdge> m_mutexSet;

    C_Progress_display my_progress_bar(
      vec_edges.size(),
      std::cout,
      "\nRelative translations computation (edge coverage algorithm)\n");

    bool bVerbose = false;

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int k = 0; k < vec_edges.size(); ++k)
    {
      #ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
      #endif
      {
        ++my_progress_bar;
      }

      const myEdge & edge = vec_edges[k];
      //-- If current edge already computed continue
      if (m_mutexSet.count(edge) || m_mutexSet.size() == vec_edges.size())
      {
        if (bVerbose)
          std::cout << "EDGES WAS PREVIOUSLY COMPUTED" << std::endl;
        continue;
      }

      std::vector<size_t> vec_possibleTriplets;
      // Find the triplets that contain the given edge
      for (size_t i = 0; i < vec_triplets.size(); ++i)
      {
        const graph::Triplet & triplet = vec_triplets[i];
        if (triplet.contain(edge))
        {
          vec_possibleTriplets.push_back(i);
        }
      }

      //-- Sort the triplet according the number of matches they have on their edges
      std::vector<size_t> vec_commonTracksPerTriplets;
      for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
      {
        vec_commonTracksPerTriplets.push_back(map_tracksPerTriplets[vec_possibleTriplets[i]]);
      }
      //-- If current edge already computed continue
      if (m_mutexSet.count(edge))
        continue;

      using namespace stl::indexed_sort;
      std::vector< sort_index_packet_descend < size_t, size_t> > packet_vec(vec_commonTracksPerTriplets.size());
      sort_index_helper(packet_vec, &vec_commonTracksPerTriplets[0]);

      std::vector<size_t> vec_possibleTripletsSorted;
      for (size_t i = 0; i < vec_commonTracksPerTriplets.size(); ++i) {
        vec_possibleTripletsSorted.push_back( vec_possibleTriplets[packet_vec[i].index] );
      }
      vec_possibleTriplets.swap(vec_possibleTripletsSorted);

      // Try to solve the triplets
      // Search the possible triplet:
      for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
      {
        const graph::Triplet & triplet = vec_triplets[vec_possibleTriplets[i]];
        //--
        // Try to estimate this triplet.
        //--
        // update precision to have good value for normalized coordinates
        double dPrecision = 8.0; // upper bound of the pixel residual
        const double ThresholdUpperBound = 1.0e-2;

        std::vector<Vec3> vec_tis(3);
        std::vector<size_t> vec_inliers;

        const std::string sOutDirectory = "./";
        const bool bTriplet_estimation = Estimate_T_triplet(
          sfm_data,
          map_globalR,
          normalized_features_provider,
          matches_provider,
          triplet,
          vec_tis,
          dPrecision,
          vec_inliers,
          ThresholdUpperBound,
          sOutDirectory);
        if (bTriplet_estimation)
        {
          std::cout << "Triplet solved: \n"
            << "#inliers: " << vec_inliers.size() << std::endl;
        }
      }
    }
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
    // Consider the pair iff it is supported by the triplet & rotation graph
    const bool b_different_pose_id = v1->id_pose != v2->id_pose;
    const int covered_pose =
      set_pose_ids.count(v1->id_pose) +
      set_pose_ids.count(v2->id_pose);
    // Different pose Id and the current edge cover the triplet edge
    if ( b_different_pose_id && covered_pose == 2 )
    {
      map_triplet_matches.insert( match_iterator );
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
  const std::vector<Mat3> vec_global_R_Triplet =
    {map_globalR.at(poses_id.i), map_globalR.at(poses_id.j), map_globalR.at(poses_id.k)};

  // check that there is enough correspondences to evaluate model
  const size_t rigSize = sfm_data.GetIntrinsics().size();
  if ( rig_tracks.size() < 50 * rigSize )
    return false ;

  // initialize rig structure for relative translation estimation
  std::vector<Vec3>  rigOffsets(rigSize);
  std::vector<Mat3>  rigRotations(rigSize);
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
  const std::map< size_t, size_t > map_poseId_to_contiguous =
    {
      {poses_id.i,0},
      {poses_id.j,1},
      {poses_id.k,2}
    };

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
      const size_t idx_view_I  = iterTrack_I->first;
      const size_t feat_I      = iterTrack_I->second;
      submapTrack::const_iterator iterTrack_J = iterTrack_I;
      std::advance(iterTrack_J, 1);

      for (; iterTrack_J != track.end(); ++iterTrack_J)
      {
        const size_t idx_view_J  = iterTrack_J->first;
        const size_t feat_J      = iterTrack_J->second;
        submapTrack::const_iterator iterTrack_K = iterTrack_J;
        std::advance(iterTrack_K, 1);

        for (; iterTrack_K != track.end(); ++iterTrack_K)
        {
          const size_t idx_view_K  = iterTrack_K->first;
          const size_t feat_K      = iterTrack_K->second;

          // initialize view structure
          const View * view_I = sfm_data.views.at(idx_view_I).get();
          const View * view_J = sfm_data.views.at(idx_view_J).get();
          const View * view_K = sfm_data.views.at(idx_view_K).get();

          // initialize intrinsic group of cameras I and J
          const IndexT intrinsic_index_I = view_I->id_intrinsic;
          const IndexT intrinsic_index_J = view_J->id_intrinsic;
          const IndexT intrinsic_index_K = view_K->id_intrinsic;

          // initialize intrinsic group of cameras I and J
          const IndexT pose_index_I = view_I->id_pose;
          const IndexT pose_index_J = view_J->id_pose;
          const IndexT pose_index_K = view_K->id_pose;
          // if the images are in 3 different poses
          if (pose_index_I != pose_index_J && pose_index_J != pose_index_K && pose_index_I != pose_index_K)
          {
            // extract normalized keypoints coordinates
            Vec2 bearing_I, bearing_J, bearing_K;
            bearing_I << normalized_features_provider->feats_per_view.at(idx_view_I).at(feat_I).coords().cast<double>();
            bearing_J << normalized_features_provider->feats_per_view.at(idx_view_J).at(feat_J).coords().cast<double>();
            bearing_K << normalized_features_provider->feats_per_view.at(idx_view_K).at(feat_K).coords().cast<double>();

            // initialize relative translation data container
            const std::vector<double> feat_cam_I = { bearing_I[0], bearing_I[1], intrinsic_index_I, map_poseId_to_contiguous.at(pose_index_I) };
            const std::vector<double> feat_cam_J = { bearing_J[0], bearing_J[1], intrinsic_index_J, map_poseId_to_contiguous.at(pose_index_J) };
            const std::vector<double> feat_cam_K = { bearing_K[0], bearing_K[1], intrinsic_index_K, map_poseId_to_contiguous.at(pose_index_K) };

            // export bearing vector in the triplet pose ordering
            std::vector<std::vector< double > > tmp(3);
            tmp[map_poseId_to_contiguous.at(pose_index_I)]= std::move(feat_cam_I);
            tmp[map_poseId_to_contiguous.at(pose_index_J)]= std::move(feat_cam_J);
            tmp[map_poseId_to_contiguous.at(pose_index_K)]= std::move(feat_cam_K);
            featsAndRigIdPerTrack.push_back( std::move(tmp) );
            sampleToTrackId[ sampleToTrackId.size() ] = cpt;
          }
        }
      }
    }
  }
  // set thresholds for relative translation estimation
  const size_t  ORSA_ITER = 1024; // max number of iterations of AC-RANSAC

  // compute model
  typedef  rigTrackTisXisTrifocalSolver  SolverType;

  typedef rig_TrackTrifocalKernel_ACRansac_N_tisXis<
    rigTrackTisXisTrifocalSolver,
    rigTrackTisXisTrifocalSolver,
    rigTrackTrifocalTensorModel> KernelType;
  KernelType kernel(featsAndRigIdPerTrack, vec_global_R_Triplet, rigRotations, rigOffsets, ThresholdUpperBound);

  rigTrackTrifocalTensorModel T;
  std::pair<double,double> acStat = robust::ACRANSAC(kernel, vec_inliers, ORSA_ITER, &T, dPrecision/minFocal, false );
  dPrecision = acStat.first;

  //-- Export data in order to have an idea of the precision of the estimates
  vec_tis = {T.t1, T.t2, T.t3};

  // update inlier list
  std::set <size_t>  inliers_tracks;
  for( size_t i = 0 ; i < vec_inliers.size() ; ++i )
      inliers_tracks.insert( sampleToTrackId[vec_inliers[i]] );

  std::vector <size_t>  inliers;
  for( size_t i = 0 ; i < rig_tracks.size() ; ++i )
     if( inliers_tracks.find(i) !=  inliers_tracks.end() )
          inliers.push_back( i );

  vec_inliers.swap( inliers );

  // if there is more than 1/3 of inliers, keep model
  const bool bTest =  ( vec_inliers.size() > 0.30 * rig_tracks.size() ) ;

  {
    std::cout << "Triplet : status: " << bTest
      << " AC: " << dPrecision
      << " inliers % " << double(inliers_tracks.size()) / rig_tracks.size() * 100.0
      << " total putative " << rig_tracks.size() << std::endl;
  }

  return bTest;

}

} // namespace sfm
} // namespace openMVG
