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
#include  "./ac_ransac_rig.hpp"
#include "./triplet_t_ACRansac_kernelAdaptator.hpp"

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/sfm/pipelines/global/mutexSet.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/linearProgramming/linearProgramming.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/sfm/pipelines/global/triplet_t_ACRansac_kernelAdaptator.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/progress/progress.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG{
namespace sfm{

// namespaces
using namespace linearProgramming;
using namespace lInfinityCV;
using namespace openMVG::trifocal;
using namespace openMVG::trifocal::kernel;
using namespace openMVG::tracks;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;


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

  // Keep the largest Biedge connected component graph of relative translations
  Pair_Set pairs;
  std::transform(
    std::begin(vec_initialRijTijEstimates), std::end(vec_initialRijTijEstimates),
    std::inserter(pairs, std::begin(pairs)),
    stl::RetrieveKey());
  const std::set<IndexT> set_remainingIds =
    openMVG::graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs, "./");
  KeepOnlyReferencedElement(set_remainingIds, vec_initialRijTijEstimates);

  std::cout << "#Remaining translations: " << vec_initialRijTijEstimates.size() << std::endl;

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
  //-------------------
  //-- GLOBAL TRANSLATIONS ESTIMATION from initial triplets t_ij guess
  //-------------------
  const std::string _sOutDirectory("./");
  {
    const std::set<IndexT> index = getIndexT(vec_initialRijTijEstimates);

    const size_t iNview = index.size();
    std::cout << "\n-------------------------------" << "\n"
      << " Global translations computation: " << "\n"
      << "   - Ready to compute " << iNview << " global translations." << "\n"
      << "     from #relative translations: " << vec_initialRijTijEstimates.size() << std::endl;

    if (iNview < 3)
    {
      // Too tiny image set to perform motion averaging
      return false;
    }
    //-- Update initial estimates from [minId,maxId] to range [0->Ncam]
    RelativeInfo_Vec vec_initialRijTijEstimates_cpy = vec_initialRijTijEstimates;
    const Pair_Set pairs = getPairs(vec_initialRijTijEstimates_cpy);
    Hash_Map<IndexT,IndexT> _reindexForward, _reindexBackward;
    reindex(pairs, _reindexForward, _reindexBackward);
    for(size_t i = 0; i < vec_initialRijTijEstimates_cpy.size(); ++i)
    {
      openMVG::relativeInfo & rel = vec_initialRijTijEstimates_cpy[i];
      rel.first = Pair(_reindexForward[rel.first.first], _reindexForward[rel.first.second]);
    }

    openMVG::system::Timer timerLP_translation;

    switch(eTranslationAveragingMethod)
    {
      case TRANSLATION_AVERAGING_L1:
      {
        double gamma = -1.0;
        std::vector<double> vec_solution;
        {
          vec_solution.resize(iNview*3 + vec_initialRijTijEstimates_cpy.size()/3 + 1);
          using namespace openMVG::linearProgramming;
          #ifdef OPENMVG_HAVE_MOSEK
            MOSEK_SolveWrapper solverLP(vec_solution.size());
          #else
            OSI_CLP_SolverWrapper solverLP(vec_solution.size());
          #endif

          lInfinityCV::Tifromtij_ConstraintBuilder_OneLambdaPerTrif cstBuilder(vec_initialRijTijEstimates_cpy);

          LP_Constraints_Sparse constraint;
          //-- Setup constraint and solver
          cstBuilder.Build(constraint);
          solverLP.setup(constraint);
          //--
          // Solving
          const bool bFeasible = solverLP.solve();
          std::cout << " \n Feasibility " << bFeasible << std::endl;
          //--
          if (bFeasible)  {
            solverLP.getSolution(vec_solution);
            gamma = vec_solution[vec_solution.size()-1];
          }
          else  {
            std::cerr << "Compute global translations: failed" << std::endl;
            return false;
          }
        }

        const double timeLP_translation = timerLP_translation.elapsed();
        //-- Export triplet statistics:
        {

          std::ostringstream os;
          os << "Translation fusion statistics.";
          os.str("");
          os << "-------------------------------" << "\n"
            << "-- #relative estimates: " << vec_initialRijTijEstimates_cpy.size()
            << " converge with gamma: " << gamma << ".\n"
            << " timing (s): " << timeLP_translation << ".\n"
            << "-------------------------------" << "\n";
          std::cout << os.str() << std::endl;
        }

        std::cout << "Found solution:\n";
        std::copy(vec_solution.begin(), vec_solution.end(), std::ostream_iterator<double>(std::cout, " "));

        std::vector<double> vec_camTranslation(iNview*3,0);
        std::copy(&vec_solution[0], &vec_solution[iNview*3], &vec_camTranslation[0]);

        std::vector<double> vec_camRelLambdas(&vec_solution[iNview*3], &vec_solution[iNview*3 + vec_initialRijTijEstimates_cpy.size()/3]);
        std::cout << "\ncam position: " << std::endl;
        std::copy(vec_camTranslation.begin(), vec_camTranslation.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\ncam Lambdas: " << std::endl;
        std::copy(vec_camRelLambdas.begin(), vec_camRelLambdas.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;

        // Update the view poses according the found camera centers
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 t(vec_camTranslation[i*3], vec_camTranslation[i*3+1], vec_camTranslation[i*3+2]);
          const IndexT camNodeId = _reindexBackward[i];
          const Mat3 & Ri = map_globalR.at(camNodeId);
          sfm_data.poses[camNodeId] = Pose3(Ri, -Ri.transpose()*t);
        }
      }
      break;

      case TRANSLATION_AVERAGING_L2:
      {
        std::vector<int> vec_edges;
        vec_edges.reserve(vec_initialRijTijEstimates_cpy.size() * 2);
        std::vector<double> vec_poses;
        vec_poses.reserve(vec_initialRijTijEstimates_cpy.size() * 3);
        std::vector<double> vec_weights;
        vec_weights.reserve(vec_initialRijTijEstimates_cpy.size());

        for(int i=0; i < vec_initialRijTijEstimates_cpy.size(); ++i)
        {
          const openMVG::relativeInfo & rel = vec_initialRijTijEstimates_cpy[i];
          vec_edges.push_back(rel.first.first);
          vec_edges.push_back(rel.first.second);
          // Since index have been remapped
          // (use the backward indexing to retrieve the second global rotation)
          const IndexT secondId = _reindexBackward[rel.first.second];
          const View * view = sfm_data.views.at(secondId).get();
          const Mat3 & Ri = map_globalR.at(view->id_pose);
          const Vec3 direction = -(Ri.transpose() * rel.second.second.normalized());

          vec_poses.push_back(direction(0));
          vec_poses.push_back(direction(1));
          vec_poses.push_back(direction(2));

          vec_weights.push_back(1.0);
        }

        const double function_tolerance = 1e-7, parameter_tolerance = 1e-8;
        const int max_iterations = 500;

        const double loss_width = 0.0; // No loss in order to compare with TRANSLATION_AVERAGING_L1

        std::vector<double> X(iNview*3, 0.0);
        if(!solve_translations_problem(
          &vec_edges[0],
          &vec_poses[0],
          &vec_weights[0],
          vec_initialRijTijEstimates_cpy.size(),
          loss_width,
          &X[0],
          function_tolerance,
          parameter_tolerance,
          max_iterations))  {
            std::cerr << "Compute global translations: failed" << std::endl;
            return false;
        }

        // Update camera center for each view
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 C(X[i*3], X[i*3+1], X[i*3+2]);
          const IndexT camNodeId = _reindexBackward[i]; // undo the reindexing
          const Mat3 & Ri = map_globalR.at(camNodeId);
          sfm_data.poses[camNodeId] = Pose3(Ri, C);
        }
      }
      break;
    }
  }
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
      //-- If current edge is already computed continue
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
        openMVG::tracks::STLMAPTracks  pose_triplet_tracks;

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
          pose_triplet_tracks,
          ThresholdUpperBound,
          sOutDirectory);

        if (bTriplet_estimation)
        {
          std::cout << "Triplet solved: \n"
            << "#inliers: " << vec_inliers.size() << std::endl;

          // Compute the 3 relative motions
          // IJ, JK, IK
          const Mat3 RI = map_globalR.at(triplet.i);
          const Mat3 RJ = map_globalR.at(triplet.j);
          const Mat3 RK = map_globalR.at(triplet.k);
          const Vec3 ti = vec_tis[0];
          const Vec3 tj = vec_tis[1];
          const Vec3 tk = vec_tis[2];

          //--- ATOMIC
          #ifdef OPENMVG_USE_OPENMP
             #pragma omp critical
          #endif
          {
            Mat3 RijGt;
            Vec3 tij;
            RelativeCameraMotion(RI, ti, RJ, tj, &RijGt, &tij);
            vec_initialEstimates.emplace_back(
              std::make_pair(triplet.i, triplet.j), std::make_pair(RijGt, tij));

            Mat3 RjkGt;
            Vec3 tjk;
            RelativeCameraMotion(RJ, tj, RK, tk, &RjkGt, &tjk);
            vec_initialEstimates.emplace_back(
              std::make_pair(triplet.j, triplet.k), std::make_pair(RjkGt, tjk));

            Mat3 RikGt;
            Vec3 tik;
            RelativeCameraMotion(RI, ti, RK, tk, &RikGt, &tik);
            vec_initialEstimates.emplace_back(
              std::make_pair(triplet.i, triplet.k), std::make_pair(RikGt, tik));

            // Add inliers as valid pairwise matches
            for (std::vector<size_t>::const_iterator iterInliers = vec_inliers.begin();
              iterInliers != vec_inliers.end(); ++iterInliers)
            {
              STLMAPTracks::const_iterator iterTracks = pose_triplet_tracks.begin();
              std::advance(iterTracks, *iterInliers);
              const submapTrack & subTrack = iterTracks->second;

              // create pairwise matches from inlier track
              for (size_t index_I = 0; index_I < subTrack.size() ; ++index_I)
              { submapTrack::const_iterator iter_I = subTrack.begin();
                std::advance(iter_I, index_I);

                // extract camera indexes
                const size_t id_view_I = iter_I->first;
                const size_t id_feat_I = iter_I->second;

                // loop on subtracks
                for (size_t index_J = index_I+1; index_J < subTrack.size() ; ++index_J)
                { submapTrack::const_iterator iter_J = subTrack.begin();
                  std::advance(iter_J, index_J);

                  // extract camera indexes
                  const size_t id_view_J = iter_J->first;
                  const size_t id_feat_J = iter_J->second;

                  newpairMatches[std::make_pair(id_view_I, id_view_J)].push_back(IndMatch(id_feat_I, id_feat_J));
                }
              }
            }
          }

          //-- Remove the 3 estimated edges
          m_mutexSet.insert(std::make_pair(triplet.i, triplet.j));
          m_mutexSet.insert(std::make_pair(triplet.j, triplet.k));
          m_mutexSet.insert(std::make_pair(triplet.i, triplet.k));
          break;
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
  openMVG::tracks::STLMAPTracks & rig_tracks,
  const double ThresholdUpperBound, //Threshold used for the trifocal tensor estimation solver used in AContrario Ransac
  const std::string & sOutDirectory) const
{
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

  // initialize rig structure for relative translation estimation
  const size_t rigSize = sfm_data.GetIntrinsics().size();
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

  // clean tracks to keep only those shared by three poses
  std::set  <size_t>          tracksToRemove;
  for (STLMAPTracks::const_iterator iterTracks = rig_tracks.begin();
    iterTracks != rig_tracks.end(); ++iterTracks)
  {
    std::set<size_t>   set_poses_index;
    // loop on subtracks
    const submapTrack & track = iterTracks->second;
    for (submapTrack::const_iterator iterTrack_I = track.begin();
        iterTrack_I != track.end(); ++iterTrack_I)
    {
      // extract pose id
      const size_t idx_view_I  = iterTrack_I->first;
      const size_t feat_I      = iterTrack_I->second;
      const View * view_I = sfm_data.views.at(idx_view_I).get();
      const IndexT pose_index_I = view_I->id_pose;

      set_poses_index.insert( pose_index_I );
    }

    // if tracks is not seen by three views, erease it
    if( set_poses_index.size() != 3 )
        tracksToRemove.insert(iterTracks->first);
  }

  // remove unneeded tracks
  for( std::set<size_t>::const_iterator iterSet = tracksToRemove.begin();
        iterSet != tracksToRemove.end(); ++iterSet)
  {
    rig_tracks.erase(*iterSet);
  }

  // check that there is enough correspondences to evaluate model
  if ( rig_tracks.size() < 50 * rigSize )
    return false ;

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
      const View * view_I = sfm_data.views.at(idx_view_I).get();
      const IndexT pose_index_I = view_I->id_pose;

      submapTrack::const_iterator iterTrack_J = iterTrack_I;
      std::advance(iterTrack_J, 1);

      for (; iterTrack_J != track.end(); ++iterTrack_J)
      {
        const size_t idx_view_J  = iterTrack_J->first;
        const size_t feat_J      = iterTrack_J->second;
        const View * view_J = sfm_data.views.at(idx_view_J).get();
        const IndexT pose_index_J = view_J->id_pose;
        if (pose_index_I == pose_index_J)
          continue;

        submapTrack::const_iterator iterTrack_K = iterTrack_J;
        std::advance(iterTrack_K, 1);

        for (; iterTrack_K != track.end(); ++iterTrack_K)
        {
          const size_t idx_view_K  = iterTrack_K->first;
          const size_t feat_K      = iterTrack_K->second;
          const View * view_K = sfm_data.views.at(idx_view_K).get();
          const IndexT pose_index_K = view_K->id_pose;

          // if the images are in 3 different poses
          if (pose_index_I != pose_index_J && pose_index_J != pose_index_K && pose_index_I != pose_index_K)
          {
            // extract normalized keypoints coordinates
            Vec2 bearing_I, bearing_J, bearing_K;
            bearing_I << normalized_features_provider->feats_per_view.at(idx_view_I).at(feat_I).coords().cast<double>();
            bearing_J << normalized_features_provider->feats_per_view.at(idx_view_J).at(feat_J).coords().cast<double>();
            bearing_K << normalized_features_provider->feats_per_view.at(idx_view_K).at(feat_K).coords().cast<double>();

            // initialize intrinsic group of cameras I and J
            const IndexT intrinsic_index_I = view_I->id_intrinsic;
            const IndexT intrinsic_index_J = view_J->id_intrinsic;
            const IndexT intrinsic_index_K = view_K->id_intrinsic;

            // initialize relative translation data container
            const std::vector<double> feat_cam_I = { bearing_I[0], bearing_I[1], intrinsic_index_I, map_poseId_to_contiguous.at(pose_index_I) };
            const std::vector<double> feat_cam_J = { bearing_J[0], bearing_J[1], intrinsic_index_J, map_poseId_to_contiguous.at(pose_index_J) };
            const std::vector<double> feat_cam_K = { bearing_K[0], bearing_K[1], intrinsic_index_K, map_poseId_to_contiguous.at(pose_index_K) };

            // export bearing vector in the triplet pose ordering
            std::vector<std::vector< double > > tmp(3);
            tmp[map_poseId_to_contiguous.at(pose_index_I)]= std::move(feat_cam_I);
            tmp[map_poseId_to_contiguous.at(pose_index_J)]= std::move(feat_cam_J);
            tmp[map_poseId_to_contiguous.at(pose_index_K)]= std::move(feat_cam_K);
            featsAndRigIdPerTrack.emplace_back( std::move(tmp) );
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
  std::pair<double,double> acStat = robust::ACRANSAC_RIG(kernel, vec_inliers, ORSA_ITER, &T, dPrecision/minFocal, false );
  // If robust estimation fail => stop.
  if (dPrecision == std::numeric_limits<double>::infinity())
    return false;

  // Update output parameters
  dPrecision = acStat.first * minFocal;

  vec_tis = {T.t1, T.t2, T.t3};

  // update inlier track list
  std::set <size_t>  inliers_tracks;
  for( size_t i = 0 ; i < vec_inliers.size() ; ++i )
      inliers_tracks.insert( sampleToTrackId[vec_inliers[i]] );

  // use move_iterator to convert the set to the vector directly ?
  std::vector <size_t>(std::make_move_iterator(inliers_tracks.begin()),
    std::make_move_iterator(inliers_tracks.end())).swap(vec_inliers);

#ifdef DEBUG_TRIPLET
  // compute 3D scene base on motion estimation
  SfM_Data    tiny_scene;

  // intialize poses (which are now shared by a group of images)
  tiny_scene.poses[poses_id.i] = Pose3(vec_global_R_Triplet[0], -vec_global_R_Triplet[0].transpose() * T.t1 );
  tiny_scene.poses[poses_id.j] = Pose3(vec_global_R_Triplet[1], -vec_global_R_Triplet[1].transpose() * T.t2 );
  tiny_scene.poses[poses_id.k] = Pose3(vec_global_R_Triplet[2], -vec_global_R_Triplet[2].transpose() * T.t3 );

  // insert views used by the relative pose pairs
  for (const auto & pairIterator : map_triplet_matches )
  {
    // initialize camera indexes
    const IndexT I = pairIterator.first.first;
    const IndexT J = pairIterator.first.second;

    // add views
    tiny_scene.views.insert(*sfm_data.GetViews().find(I));
    tiny_scene.views.insert(*sfm_data.GetViews().find(J));

    // add intrinsics
    const View * view_I = sfm_data.GetViews().at(I).get();
    const View * view_J = sfm_data.GetViews().at(J).get();
    tiny_scene.intrinsics.insert(*sfm_data.GetIntrinsics().find(view_I->id_intrinsic));
    tiny_scene.intrinsics.insert(*sfm_data.GetIntrinsics().find(view_J->id_intrinsic));
  }

  // Fill sfm_data with the inliers tracks. Feed image observations: no 3D yet.
  Landmarks & structure = tiny_scene.structure;
  for (size_t idx=0; idx < vec_inliers.size(); ++idx)
  {
    const size_t trackId = vec_inliers[idx];
    const submapTrack & track = rig_tracks[trackId];
    Observations & obs = structure[idx].obs;
    for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
    {
      // get view Id and feat ID
      const size_t viewIndex = it->first;
      const size_t featIndex = it->second;

      // initialize view and get intrinsics
      const View * view = sfm_data.GetViews().at(viewIndex).get();
      const cameras::IntrinsicBase *  cam = sfm_data.GetIntrinsics().find(view->id_intrinsic)->second.get();
      const cameras::Rig_Pinhole_Intrinsic * rig_intrinsicPtr = dynamic_cast< const cameras::Rig_Pinhole_Intrinsic * >(cam);
      const Vec2  principal_point = rig_intrinsicPtr->principal_point();

      // get normalized feature
      const PointFeature & pt = normalized_features_provider->feats_per_view.at(viewIndex)[featIndex];
      PointFeature pt_unnormalized( pt.x() * rig_intrinsicPtr->focal() + principal_point.x(),
                                    pt.y() * rig_intrinsicPtr->focal() + principal_point.y());

      obs[viewIndex] = Observation(pt_unnormalized.coords().cast<double>(), featIndex);
    }
  }

  // Compute 3D landmark positions (triangulation of the observations)
  {
    SfM_Data_Structure_Computation_Blind structure_estimator(false);
    structure_estimator.triangulate(tiny_scene);
  }

  // export scene for visualization
  std::ostringstream os;
  os << poses_id.i << "_" << poses_id.j << "_" << poses_id.k << ".ply";
  Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
#endif

  // if there is more than 1/3 of inliers, keep model
  const bool bTest =  ( vec_inliers.size() > 0.66 * rig_tracks.size() ) ;

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
