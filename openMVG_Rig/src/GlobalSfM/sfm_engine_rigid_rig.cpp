/*
 * sfm_engine_rigid_rig
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

#include "./sfm_engine_rigid_rig.hpp"
#include "./sfm_engine_translation_averaging.hpp"
#include "./sfm_robust_relative_pose_rig.hpp"
#include "openMVG/multiview/essential.hpp"

#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/stl/stl.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;

ReconstructionEngine_RelativeMotions_RigidRig::ReconstructionEngine_RelativeMotions_RigidRig(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory), _sLoggingFile(sloggingFile), _normalized_features_provider(NULL) {

  if (!_sLoggingFile.empty())
  {
    // setup HTML logger
    _htmlDocStream = std::make_shared<htmlDocument::htmlDocumentStream>("GlobalReconstructionEngine SFM report.");
    _htmlDocStream->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("ReconstructionEngine_RelativeMotions_RigidRig")));
    _htmlDocStream->pushInfo("<hr>");

    _htmlDocStream->pushInfo( "Dataset info:");
    _htmlDocStream->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }

  // Set default motion Averaging methods
  _eRotationAveragingMethod = ROTATION_AVERAGING_L2;
  _eTranslationAveragingMethod = TRANSLATION_AVERAGING_L1;
}

ReconstructionEngine_RelativeMotions_RigidRig::~ReconstructionEngine_RelativeMotions_RigidRig()
{
  if (!_sLoggingFile.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(_sLoggingFile.c_str());
    htmlFileStream << _htmlDocStream->getDoc();
  }
}

void ReconstructionEngine_RelativeMotions_RigidRig::SetFeaturesProvider(Features_Provider * provider)
{
  _features_provider = provider;

  // Copy features and save a normalized version
  _normalized_features_provider = std::make_shared<Features_Provider>(*provider);
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (Hash_Map<IndexT, PointFeatures>::iterator iter = _normalized_features_provider->feats_per_view.begin();
    iter != _normalized_features_provider->feats_per_view.end(); ++iter)
  {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
#endif
    {
      // get the related view & camera intrinsic and compute the corresponding bearing vectors
      const View * view = _sfm_data.GetViews().at(iter->first).get();
      const std::shared_ptr<IntrinsicBase> cam = _sfm_data.GetIntrinsics().find(view->id_intrinsic)->second;
      for (PointFeatures::iterator iterPt = iter->second.begin();
        iterPt != iter->second.end(); ++iterPt)
      {
        const Vec3 bearingVector = (*cam)(cam->get_ud_pixel(iterPt->coords().cast<double>()));
        iterPt->coords() << (bearingVector.head(2) / bearingVector(2)).cast<float>();
      }
    }
  }
}

void ReconstructionEngine_RelativeMotions_RigidRig::SetMatchesProvider(Matches_Provider * provider)
{
  _matches_provider = provider;
}

void ReconstructionEngine_RelativeMotions_RigidRig::SetRotationAveragingMethod
(
  ERotationAveragingMethod eRotationAveragingMethod
)
{
  _eRotationAveragingMethod = eRotationAveragingMethod;
}

void ReconstructionEngine_RelativeMotions_RigidRig::SetTranslationAveragingMethod
(
  ETranslationAveragingMethod eTranslationAveragingMethod
)
{
  _eTranslationAveragingMethod = eTranslationAveragingMethod;
}

bool ReconstructionEngine_RelativeMotions_RigidRig::Process() {

  //-------------------
  // TODO: keep largest biedge connected pose subgraph
  //-------------------

  openMVG::rotation_averaging::RelativeRotations relatives_R;
  Compute_Relative_Rotations(relatives_R);

  Hash_Map<IndexT, Mat3> global_rotations;
  if (!Compute_Global_Rotations(relatives_R, global_rotations))
  {
    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
    return false;
  }
  matching::PairWiseMatches  tripletWise_matches;
  if (!Compute_Global_Translations(global_rotations, relatives_R, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!_sLoggingFile.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    _htmlDocStream->pushInfo("<hr>");
    _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << _sfm_data.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << _sfm_data.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << _sfm_data.GetPoses().size() << "<br>"
      << "-- Track count: "  << _sfm_data.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    _htmlDocStream->pushInfo(os.str());
  }

  return true;
}

/// Compute from relative rotations the global rotations of the camera poses
bool ReconstructionEngine_RelativeMotions_RigidRig::Compute_Global_Rotations
(
  const rotation_averaging::RelativeRotations & relatives_R,
  Hash_Map<IndexT, Mat3> & global_rotations
)
{
  if(relatives_R.empty())
    return false;
  // Log statistics about the relative rotation graph
  {
    std::set<IndexT> set_pose_ids;
    for (const auto & relative_R : relatives_R)
    {
      set_pose_ids.insert(relative_R.i);
      set_pose_ids.insert(relative_R.j);
    }

    std::cout << "\n-------------------------------" << "\n"
      << " Global rotations computation: " << "\n"
      << "  #relative rotations: " << relatives_R.size() << "\n"
      << "  #global rotations: " << set_pose_ids.size() << std::endl;
  }

  // Global Rotation solver:
  const ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod =
    TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR;
    //TRIPLET_ROTATION_INFERENCE_NONE;

  GlobalSfM_Rotation_AveragingSolver rotation_averaging_solver;
  const bool b_rotation_averaging = rotation_averaging_solver.Run(
    _eRotationAveragingMethod, eRelativeRotationInferenceMethod,
    relatives_R, global_rotations);

  std::cout << "Found #global_rotations: " << global_rotations.size() << std::endl;

  if (b_rotation_averaging)
  {
    // Log input graph to the HTML report
    if (!_sLoggingFile.empty() && !_sOutDirectory.empty())
    {
      // Log a relative pose graph
      {
        std::set<IndexT> set_pose_ids;
        Pair_Set relative_pose_pairs;
        for (const auto & view : _sfm_data.GetViews())
        {
          const IndexT pose_id = view.second->id_pose;
          set_pose_ids.insert(pose_id);
        }
        const std::string sGraph_name = "global_relative_rotation_pose_graph_final";
        graph::indexedGraph putativeGraph(set_pose_ids, rotation_averaging_solver.GetUsedPairs());
        graph::exportToGraphvizData(
          stlplus::create_filespec(_sOutDirectory, sGraph_name),
          putativeGraph.g);

        using namespace htmlDocument;
        std::ostringstream os;

        os << "<br>" << sGraph_name << "<br>"
           << "<img src=\""
           << stlplus::create_filespec(_sOutDirectory, sGraph_name, "svg")
           << "\" height=\"600\">\n";

        _htmlDocStream->pushInfo(os.str());
      }
    }
  }

  return b_rotation_averaging;
}

/// Compute/refine relative translations and compute global translations
bool ReconstructionEngine_RelativeMotions_RigidRig::Compute_Global_Translations
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  const openMVG::rotation_averaging::RelativeRotations & relatives_R,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Translation averaging (compute translations & update them to a global common coordinates system)
  GlobalSfMRig_Translation_AveragingSolver translation_averaging_solver;
  const bool bTranslationAveraging = translation_averaging_solver.Run(
    _eTranslationAveragingMethod,
    _sfm_data,
    _normalized_features_provider.get(),
    _matches_provider,
    global_rotations,
    relatives_R,
    tripletWise_matches);

  if (!_sLoggingFile.empty())
  {
    Save(_sfm_data,
      stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return bTranslationAveraging;
}

/// Compute the initial structure of the scene
bool ReconstructionEngine_RelativeMotions_RigidRig::Compute_Initial_Structure
(
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Build tracks from selected triplets (Union of all the validated triplet tracks (_tripletWise_matches))
  {
    using namespace openMVG::tracks;
    TracksBuilder tracksBuilder;
    tracksBuilder.Build(tripletWise_matches);
    tracksBuilder.Filter(3);
    STLMAPTracks map_selectedTracks; // reconstructed track (visibility per 3D point)
    tracksBuilder.ExportToSTL(map_selectedTracks);

    // Fill sfm_data with the computed tracks (no 3D yet)
    Landmarks & structure = _sfm_data.structure;
    IndexT idx(0);
    for (STLMAPTracks::const_iterator itTracks = map_selectedTracks.begin();
      itTracks != map_selectedTracks.end();
      ++itTracks, ++idx)
    {
      const submapTrack & track = itTracks->second;
      structure[idx] = Landmark();
      Observations & obs = structure.at(idx).obs;
      for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
      {
        const size_t imaIndex = it->first;
        const size_t featIndex = it->second;
        const PointFeature & pt = _features_provider->feats_per_view.at(imaIndex)[featIndex];
        obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
      }
    }

    std::cout << std::endl << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats:
      //    - number of images
      //    - number of tracks
      std::set<size_t> set_imagesId;
      TracksUtilsMap::ImageIdInTracks(map_selectedTracks, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<size_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<size_t, size_t> map_Occurence_TrackLength;
      TracksUtilsMap::TracksLength(map_selectedTracks, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (std::map<size_t, size_t>::const_iterator iter = map_Occurence_TrackLength.begin();
        iter != map_Occurence_TrackLength.end(); ++iter)  {
        osTrack << "\t" << iter->first << "\t" << iter->second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }

  // Compute 3D position of the landmark of the structure by triangulation of the observations
  {
    openMVG::system::Timer timer;

    const IndexT trackCountBefore = _sfm_data.GetLandmarks().size();
    SfM_Data_Structure_Computation_Blind structure_estimator(true);
    structure_estimator.triangulate(_sfm_data);

    std::cout << "\n#removed tracks (invalid triangulation): " <<
      trackCountBefore - IndexT(_sfm_data.GetLandmarks().size()) << std::endl;
    std::cout << std::endl << "  Triangulation took (s): " << timer.elapsed() << std::endl;

    // Export initial structure
    if (!_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "initial_structure", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  return !_sfm_data.structure.empty();
}

// Adjust the scene (& remove outliers)
bool ReconstructionEngine_RelativeMotions_RigidRig::Adjust()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, false, true, false);
  if (b_BA_Status)
  {
    if (!_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, true, true, false);
    if (b_BA_Status && !_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && !_bFixedIntrinsics) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, true, true, true);
    if (b_BA_Status && !_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = _sfm_data.structure.size();
  RemoveOutliers_PixelResidualError(_sfm_data, 4.0);
  const size_t pointcount_pixelresidual_filter = _sfm_data.structure.size();
  RemoveOutliers_AngleError(_sfm_data, 2.0);
  const size_t pointcount_angular_filter = _sfm_data.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!_sLoggingFile.empty())
  {
    Save(_sfm_data,
      stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, true, true, false);
  if (b_BA_Status && !_sLoggingFile.empty())
  {
    Save(_sfm_data,
      stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_04_refined_RT_Xi", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLenght = 3; // 2 min
  if (eraseUnstablePosesAndObservations(_sfm_data, minPointPerPose, minTrackLenght))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = _sfm_data.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";

    b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, true, true, !_bFixedIntrinsics);
    if (b_BA_Status && !_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_05_unstable_removed", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  return b_BA_Status;
}

void ReconstructionEngine_RelativeMotions_RigidRig::Compute_Relative_Rotations
(
  rotation_averaging::RelativeRotations & vec_relatives_R
)
{
  //
  // Build the Relative pose graph from matches:
  //
  /// pairwise view relation shared between poseIds
  typedef std::map< Pair, Pair_Set > RigWiseMatches;

  // List shared correspondences (pairs) between poses
  RigWiseMatches rigWiseMatches;
  for (PairWiseMatches::const_iterator iterMatches = _matches_provider->_pairWise_matches.begin();
    iterMatches != _matches_provider->_pairWise_matches.end(); ++iterMatches)
  {
    const Pair pair = iterMatches->first;
    const View * v1 = _sfm_data.GetViews().at(pair.first).get();
    const View * v2 = _sfm_data.GetViews().at(pair.second).get();
    rigWiseMatches[Pair(v1->id_pose, v2->id_pose)].insert(pair);
  }

  C_Progress_display my_progress_bar( rigWiseMatches.size(),
      std::cout, "\n- Relative pose computation -\n" );

#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  // For each non-central camera pairs, compute the relative pose from pairwise point matches:
  for (const auto & relativePoseIterator : rigWiseMatches)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      const Pair relative_pose_pair = relativePoseIterator.first;
      const Pair_Set & match_pairs = relativePoseIterator.second;

      // If a pair has the same ID, discard it
      if ( relative_pose_pair.first != relative_pose_pair.second )
      {
        // Compute the relative pose...
        // - if many pairs: relative pose between rigid rigs
        // TODO -> if only one pair of matches: relative pose between two pinhole images

        //--
        // Setup intrinsics data used by the relative poses
        // - focal
        // - subposes (if any)
        //--

        // create rig structure using openGV
        opengv::translations_t  rigOffsets;
        opengv::rotations_t     rigRotations;
        double                  minFocal = 1.0e10;
        std::map<IndexT, IndexT> intrinsic_id_remapping; // to ensure contiguous indices

        // Update rig structure from OpenMVG data to OpenGV convention
        // - List all used views, for each parse and use intrinsic data
        std::set<IndexT> used_views;
        IndexT intrinsic_number_pose_one=0;
        IndexT intrinsic_number_pose_two=0;
        for (const auto & pairIterator : match_pairs )
        {
          const IndexT I = pairIterator.first;
          const IndexT J = pairIterator.second;

          const View * view_I = _sfm_data.views[I].get();
          const View * view_J = _sfm_data.views[J].get();

          // Check that valid camera are existing for the pair view
          if (_sfm_data.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
            _sfm_data.GetIntrinsics().count(view_J->id_intrinsic) == 0)
            continue;

          if (intrinsic_id_remapping.count(view_I->id_intrinsic) == 0)
          {
            used_views.insert(I);
            intrinsic_id_remapping.insert(std::make_pair(view_I->id_intrinsic, intrinsic_id_remapping.size()));

            // update number of subcameras for each pose
            if(view_I->id_pose == relative_pose_pair.first )
              ++intrinsic_number_pose_one;
            else
              ++intrinsic_number_pose_two;
          }
          if (intrinsic_id_remapping.count(view_J->id_intrinsic) == 0)
          {
            used_views.insert(J);
            intrinsic_id_remapping.insert(std::make_pair(view_J->id_intrinsic, intrinsic_id_remapping.size()));

            // update number of subcamera for each poses
            if(view_J->id_pose == relative_pose_pair.first )
              ++intrinsic_number_pose_one;
            else
              ++intrinsic_number_pose_two;
          }
        }
        rigOffsets.resize(intrinsic_id_remapping.size());
        rigRotations.resize(intrinsic_id_remapping.size());
        for (const auto & view_id : used_views)
        {
          const View * view = _sfm_data.views[view_id].get();
          const cameras::IntrinsicBase * intrinsicPtr = _sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
          if ( intrinsicPtr->getType() == cameras::PINHOLE_RIG_CAMERA )
          {
            // retrieve camera information from shared pointer
            const cameras::Rig_Pinhole_Intrinsic * rig_intrinsicPtr = dynamic_cast< const cameras::Rig_Pinhole_Intrinsic * > (intrinsicPtr);
            const geometry::Pose3 sub_pose = rig_intrinsicPtr->get_subpose();

            // update rig stucture
            const IndexT index  = intrinsic_id_remapping.at(view->id_intrinsic);
            rigOffsets[index]   = sub_pose.center();
            rigRotations[index] = sub_pose.rotation().transpose();

            minFocal = std::min( minFocal , rig_intrinsicPtr->focal() );
          }
        }

        //--
        // Setup bearing vector correspondences by using feature tracking
        //--

        // Select the matches that are shared by view having the selected pose Ids
        // Will be used to compute tracks
        matching::PairWiseMatches  matches_rigpair;
        for (const auto & pairIterator : match_pairs )
        {
          const IndexT I = pairIterator.first;
          const IndexT J = pairIterator.second;

          const View * view_I = _sfm_data.views[I].get();
          const View * view_J = _sfm_data.views[J].get();

          // Check that valid camera are existing for the pair view
          if (_sfm_data.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
            _sfm_data.GetIntrinsics().count(view_J->id_intrinsic) == 0)
            continue;

          // add bearing vectors if they do not belong to the same pose
          if ( view_I->id_pose == view_J->id_pose )
            continue;

          // export pairwise matches
          matches_rigpair.insert( std::make_pair( pairIterator, _matches_provider->_pairWise_matches.at( pairIterator )));
        }

        // initialize tracks
        using namespace openMVG::tracks;
        TracksBuilder tracksBuilder;
        tracksBuilder.Build( matches_rigpair );
        tracksBuilder.Filter( 2 ); // matches must be seen by at least two view/pose.
        STLMAPTracks map_tracks; // reconstructed track (visibility per 3D point)
        tracksBuilder.ExportToSTL(map_tracks);

        if ( map_tracks.size() > 50 * std::max(intrinsic_number_pose_one, intrinsic_number_pose_two)) // Early rejection if too few tracks
        {
          // initialize structures used for matching between rigs
          opengv::bearingVectors_t bearingVectorsRigOne, bearingVectorsRigTwo;
          std::vector<int> camCorrespondencesRigOne, camCorrespondencesRigTwo;
          bearingVectorsRigOne.reserve(map_tracks.size());
          bearingVectorsRigTwo.reserve(map_tracks.size());
          camCorrespondencesRigOne.reserve(map_tracks.size());
          camCorrespondencesRigTwo.reserve(map_tracks.size());

          // List the tracks to associate a pair of bearing vector to a track Id
          std::vector<size_t>  vec_bearingIdToTrackId;
          vec_bearingIdToTrackId.reserve(map_tracks.size());
          for (STLMAPTracks::const_iterator iterTracks = map_tracks.begin();
            iterTracks != map_tracks.end(); ++iterTracks)
          {
            const submapTrack & track = iterTracks->second;
            for (submapTrack::const_iterator iterTrack_I = track.begin();
              iterTrack_I != track.end(); ++iterTrack_I)
            {
              const size_t I  = iterTrack_I->first;
              const size_t feat_I = iterTrack_I->second;
              submapTrack::const_iterator iterTrack_J = iterTrack_I;
              std::advance(iterTrack_J, 1);
              for (; iterTrack_J != track.end(); ++iterTrack_J)
              {
                const size_t J  = iterTrack_J->first;
                const size_t feat_J = iterTrack_J->second;

                // initialize view structure
                const View * view_I = _sfm_data.views[I].get();
                const View * view_J = _sfm_data.views[J].get();

                if (view_I->id_pose == view_J->id_pose)
                {
                  continue;
                  // This correspondences is meaningless for relative pose estimation,
                  // it belongs to the same rig.
                }

                // Add bearing_vector correspondences to the list
                opengv::bearingVector_t  bearing_vector;

                // extract normalized keypoints coordinates
                bearing_vector << _normalized_features_provider->feats_per_view[I][feat_I].coords().cast<double>(), 1.0;
                bearing_vector.normalize();
                bearingVectorsRigOne.push_back( bearing_vector );
                camCorrespondencesRigOne.push_back( intrinsic_id_remapping.at(view_I->id_intrinsic) );

                bearing_vector << _normalized_features_provider->feats_per_view[J][feat_J].coords().cast<double>(), 1.0;
                bearing_vector.normalize();
                bearingVectorsRigTwo.push_back( bearing_vector );
                camCorrespondencesRigTwo.push_back( intrinsic_id_remapping.at(view_J->id_intrinsic) );

                // update map
                vec_bearingIdToTrackId.push_back(iterTracks->first);
              }
            }
          }// end loop on tracks

          //--
          // Robust pose estimation
          //--

          //--> Estimate the best possible Rotation/Translation from correspondences
          double errorMax = std::numeric_limits<double>::max();
          const double upper_bound_pixel_threshold = 2.5 * sqrt(2.0);
          const double maxExpectedError = 1.0 - cos ( atan ( upper_bound_pixel_threshold / minFocal ) );

          opengv::transformation_t pose;
          std::vector<size_t> vec_inliers;

          const bool isPoseUsable = SfMRobust::robustRigPose(
            bearingVectorsRigOne, bearingVectorsRigTwo,
            camCorrespondencesRigOne, camCorrespondencesRigTwo,
            rigOffsets, rigRotations,
            &pose,
            &vec_inliers,
            &errorMax,
            maxExpectedError);

          if ( isPoseUsable )
          {
            // set output model
            geometry::Pose3 relativePose(pose.block<3,3>(0,0).transpose(), pose.col(3));

            // Build a tiny SfM scene with only the geometry of the relative pose
            //  for parameters refinement: 3D points & camera poses.
            SfM_Data tiny_scene;

            const IndexT indexRig1 = relative_pose_pair.first;
            const IndexT indexRig2 = relative_pose_pair.second;
            // intialize poses (shared by a group of images)
            tiny_scene.poses[indexRig1] = Pose3(Mat3::Identity(), Vec3::Zero());
            tiny_scene.poses[indexRig2] = relativePose;

            // insert views used by the relative pose pairs
            for (const auto & pairIterator : match_pairs)
            {
              // initialize camera indexes
              const IndexT I = pairIterator.first;
              const IndexT J = pairIterator.second;

              // add views
              tiny_scene.views.insert(*_sfm_data.GetViews().find(pairIterator.first));
              tiny_scene.views.insert(*_sfm_data.GetViews().find(pairIterator.second));

              // add intrinsics
              const View * view_I = _sfm_data.views[I].get();
              const View * view_J = _sfm_data.views[J].get();
              tiny_scene.intrinsics.insert(*_sfm_data.GetIntrinsics().find(view_I->id_intrinsic));
              tiny_scene.intrinsics.insert(*_sfm_data.GetIntrinsics().find(view_J->id_intrinsic));
            }

            // Fill sfm_data with the inliers tracks. Feed image observations: no 3D yet.
            Landmarks & structure = tiny_scene.structure;
            for (size_t idx=0; idx < vec_inliers.size(); ++idx)
            {
              const size_t trackId = vec_bearingIdToTrackId[vec_inliers[idx]];
              const submapTrack & track = map_tracks.at(trackId);
              Observations & obs = structure[idx].obs;
              for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
              {
                const size_t imaIndex = it->first;
                const size_t featIndex = it->second;
                const PointFeature & pt = _features_provider->feats_per_view.at(imaIndex)[featIndex];
                obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
              }
            }

            // Compute 3D position of the landmarks (triangulation of the observations)
            {
              SfM_Data_Structure_Computation_Blind structure_estimator(false);
              structure_estimator.triangulate(tiny_scene);
            }

            // Refine structure and poses (keep intrinsic constant)
            Bundle_Adjustment_Ceres::BA_options options(false, false);
            options._linear_solver_type = ceres::DENSE_SCHUR;
            Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
            if (bundle_adjustment_obj.Adjust(tiny_scene, true, true, false))
            {
              // --> to debug: save relative pair geometry on disk
              //std::ostringstream os;
              //os << relative_pose_pair.first << "_" << relative_pose_pair.second << ".ply";
              //Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
              //
              const Mat3 R1 = tiny_scene.poses[indexRig1].rotation();
              const Mat3 R2 = tiny_scene.poses[indexRig2].rotation();
              const Vec3 t1 = tiny_scene.poses[indexRig1].translation();
              const Vec3 t2 = tiny_scene.poses[indexRig2].translation();
              // Compute relative motion and save it
              Mat3 Rrel;
              Vec3 trel;
              RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
              // Update found relative pose
              relativePose = Pose3(Rrel, -Rrel.transpose() * trel);

              // Add the relative rotation to the relative 'rotation' pose graph
              using namespace openMVG::rotation_averaging;
              vec_relatives_R.push_back(RelativeRotation(
                relative_pose_pair.first, relative_pose_pair.second,
                relativePose.rotation(), vec_inliers.size() ));
            }
          }
        }
      }
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
      {
        ++my_progress_bar;
      }
    } // omp single no wait
  }
  // Log input graph to the HTML report
  if (!_sLoggingFile.empty() && !_sOutDirectory.empty())
  {
    // Log a relative view graph
    {
      std::set<IndexT> set_ViewIds;
      std::transform(_sfm_data.GetViews().begin(), _sfm_data.GetViews().end(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(_matches_provider->_pairWise_matches));
      graph::exportToGraphvizData(
        stlplus::create_filespec(_sOutDirectory, "global_relative_rotation_view_graph"),
        putativeGraph.g);
    }

    // Log a relative pose graph
    {
      std::set<IndexT> set_pose_ids;
      Pair_Set relative_pose_pairs;
      for (const auto & relative_R : vec_relatives_R)
      {
        const Pair relative_pose_indices(relative_R.i, relative_R.j);
        relative_pose_pairs.insert(relative_pose_indices);
        set_pose_ids.insert(relative_R.i);
        set_pose_ids.insert(relative_R.j);
      }
      const std::string sGraph_name = "global_relative_rotation_pose_graph";
      graph::indexedGraph putativeGraph(set_pose_ids, relative_pose_pairs);
      graph::exportToGraphvizData(
        stlplus::create_filespec(_sOutDirectory, sGraph_name),
        putativeGraph.g);
            using namespace htmlDocument;
      std::ostringstream os;

      os << "<br>" << "global_relative_rotation_pose_graph" << "<br>"
         << "<img src=\""
         << stlplus::create_filespec(_sOutDirectory, "global_relative_rotation_pose_graph", "svg")
         << "\" height=\"600\">\n";

      _htmlDocStream->pushInfo(os.str());
    }
  }
}

} // namespace sfm
} // namespace openMVG
