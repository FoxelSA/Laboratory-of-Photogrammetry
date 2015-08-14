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
#include "./sfm_robust_relative_pose_rig.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/multiview/essential.hpp"

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
  for (Hash_Map<IndexT, PointFeatures>::iterator iter = _normalized_features_provider->feats_per_view.begin();
    iter != _normalized_features_provider->feats_per_view.end(); ++iter)
  {
    // get the related view & camera intrinsic and compute the corresponding bearing vectors
    const View * view = _sfm_data.GetViews().at(iter->first).get();
    const std::shared_ptr<IntrinsicBase> cam = _sfm_data.GetIntrinsics().find(view->id_intrinsic)->second;
    for (PointFeatures::iterator iterPt = iter->second.begin();
      iterPt != iter->second.end(); ++iterPt)
    {
      const Vec3 bearingVector = (*cam)(cam->get_ud_pixel(iterPt->coords().cast<double>()));
      const Vec2 bearingVectorNormalized = bearingVector.head(2) / bearingVector(2);
      iterPt->coords() = Vec2f(bearingVectorNormalized(0), bearingVectorNormalized(1));
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

  RelativeInfo_Map relative_Rt;
  Compute_Relative_Rotations(relative_Rt);

  if (!Compute_Global_Rotations())
  {
    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Global_Translations())
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Initial_Structure())
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
bool ReconstructionEngine_RelativeMotions_RigidRig::Compute_Global_Rotations()
{
  return false;
}

/// Compute/refine relative translations and compute global translations
bool ReconstructionEngine_RelativeMotions_RigidRig::Compute_Global_Translations()
{
  return false;
}

/// Compute the initial structure of the scene
bool ReconstructionEngine_RelativeMotions_RigidRig::Compute_Initial_Structure()
{
  return false;
}

// Adjust the scene (& remove outliers)
bool ReconstructionEngine_RelativeMotions_RigidRig::Adjust()
{
  return false;
}

void ReconstructionEngine_RelativeMotions_RigidRig::Compute_Relative_Rotations(RelativeInfo_Map & vec_relatives)
{
  //
  // Build the Relative pose graph:
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

  // create rig structure using openGV
  opengv::translations_t  rigOffsets;
  opengv::rotations_t     rigRotations;
  double                  minFocal=1.0e10;

  //update size
  rigOffsets.resize  ( _sfm_data.GetIntrinsics().size() );
  rigRotations.resize( _sfm_data.GetIntrinsics().size() );

  // affect rig structure using opengv
  for (const auto & intrinsicVal : _sfm_data.GetIntrinsics())
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
      rigRotations[index] = sub_pose.rotation().transpose();

      minFocal = std::min( minFocal , focal );
    }
  }

  // For each non-central camera pairs, compute the rotation from pairwise point matches:
  for (const auto & relativePoseIterator : rigWiseMatches)
  {
    const Pair relative_pose_pair = relativePoseIterator.first;
    const Pair_Set & match_pairs = relativePoseIterator.second;

    const IndexT   indexRig1 = relative_pose_pair.first;
    const IndexT   indexRig2 = relative_pose_pair.second;

    //if the consider pair has the same ID, discard it
    if ( relative_pose_pair.first != relative_pose_pair.second )
    {
      // copy to a vector for use with threading
      const Pair_Vec pair_vec(match_pairs.begin(), match_pairs.end());

      // Compute the relative pose...
      // - if only one pair of matches: relative pose between two pinhole images
      // - if many pairs: relative pose between rigid rigs

      // initialize structure used for matching between rigs
      opengv::bearingVectors_t bearingVectorsRigOne;
      opengv::bearingVectors_t bearingVectorsRigTwo;

      std::vector<int>  camCorrespondencesRigOne;
      std::vector<int>  camCorrespondencesRigTwo;

      opengv::transformation_t  pose;
      std::vector<size_t> vec_inliers;

      // initialise matches between the two rigs
      matching::PairWiseMatches  matches_rigpair;

      // create pairwise matches in order to compute tracks
      for (const auto & pairIterator : match_pairs )
      {
        const size_t I = pairIterator.first;
        const size_t J = pairIterator.second;

        const View * view_I = _sfm_data.views[I].get();
        const View * view_J = _sfm_data.views[J].get();

        // Check that valid camera are existing for the pair index
        if (_sfm_data.GetIntrinsics().find(view_I->id_intrinsic) == _sfm_data.GetIntrinsics().end() ||
          _sfm_data.GetIntrinsics().find(view_J->id_intrinsic) == _sfm_data.GetIntrinsics().end())
          continue;

        // export pairwise matches
        matches_rigpair[pairIterator] = _matches_provider->_pairWise_matches.at( pairIterator );
      }

      // initialize tracks
      using namespace openMVG::tracks;
      TracksBuilder tracksBuilder;
      tracksBuilder.Build( matches_rigpair );
      tracksBuilder.Filter( 2 );
      STLMAPTracks map_tracks; // reconstructed track (visibility per 3D point)
      tracksBuilder.ExportToSTL(map_tracks);

      // extract associated subcamera id for each tracks
      std::map < size_t, size_t >  map_bearingIdToTrackId;
      size_t cpt = 0;
      for (STLMAPTracks::const_iterator iterTracks = map_tracks.begin();
        iterTracks != map_tracks.end(); ++iterTracks, ++cpt)
      {
        const submapTrack & subTrack = iterTracks->second;
        std::vector<std::pair<size_t, size_t> >  imgAndFeat_rigI;
        std::vector<std::pair<size_t, size_t> >  imgAndFeat_rigJ;

        for (submapTrack::const_iterator iterSubTrack = subTrack.begin();
                iterSubTrack != subTrack.end(); ++iterSubTrack)
        {
          const size_t imaIndex  = iterSubTrack->first;
          const size_t featIndex = iterSubTrack->second;

          const View * view = _sfm_data.views[imaIndex].get();
          const size_t rigidId = view->id_pose;

          if( rigidId == indexRig1 )
              imgAndFeat_rigI.push_back(std::make_pair(imaIndex,featIndex));

          if( rigidId == indexRig2 )
              imgAndFeat_rigJ.push_back(std::make_pair(imaIndex,featIndex));
        }

        // loop on matches between rigs
        for( size_t  indI = 0; indI < imgAndFeat_rigI.size(); ++indI )
        {
            for( size_t indJ = 0; indJ < imgAndFeat_rigJ.size(); ++indJ )
            {
              // extract image index and feat index
              const size_t  I      = imgAndFeat_rigI[indI].first;
              const size_t  feat_I = imgAndFeat_rigI[indI].second;

              const size_t  J      = imgAndFeat_rigJ[indJ].first;
              const size_t  feat_J = imgAndFeat_rigJ[indJ].second;

              // initialize view structure
              const View * view_I = _sfm_data.views[I].get();
              const View * view_J = _sfm_data.views[J].get();

              // initialise intrinsic group of cameras I and J
              const IndexT intrinsic_index_I = view_I->id_intrinsic;
              const IndexT intrinsic_index_J = view_J->id_intrinsic;

              // extract features
              opengv::bearingVector_t  bearing_one;
              opengv::bearingVector_t  bearing_two;

              // extract normalized keypoints coordinates
              bearing_one << _normalized_features_provider->feats_per_view[I][feat_I].coords().cast<double>(), 1.0;
              bearing_two << _normalized_features_provider->feats_per_view[J][feat_J].coords().cast<double>(), 1.0;

              // normalize bearing vectors
              bearing_one.normalized();
              bearing_two.normalized();

              // add bearing vectors to list and update correspondences list
              bearingVectorsRigOne.push_back( bearing_one );
              camCorrespondencesRigOne.push_back( intrinsic_index_I );

              // add bearing vectors to list and update correspondences list
              bearingVectorsRigTwo.push_back( bearing_two );
              camCorrespondencesRigTwo.push_back( intrinsic_index_J );

              // update map
              map_bearingIdToTrackId[bearingVectorsRigTwo.size()-1] = cpt;
            }
        }

      }// end loop on tracks

      //--> Estimate the best possible Rotation/Translation from correspondences
      double errorMax = std::numeric_limits<double>::max();
      const double maxExpectedError = 1.0 - cos ( atan ( sqrt(2.0) * 2.5 / minFocal ) );
      bool isPoseUsable = false ;

      if ( map_tracks.size() > 50 * rigOffsets.size())
      {
          isPoseUsable = SfMRobust::robustRigPose(
                                bearingVectorsRigOne,
                                bearingVectorsRigTwo,
                                camCorrespondencesRigOne,
                                camCorrespondencesRigTwo,
                                rigOffsets,
                                rigRotations,
                                &pose,
                                &vec_inliers,
                                &errorMax,
                                maxExpectedError);
      }

      if( isPoseUsable )
      {
        // set output model
        geometry::Pose3  relativePose( pose.block<3,3>(0,0).transpose(), pose.col(3));

        // keep only tracks related to inliers
        std::set<size_t>   inliers_set;
        for(size_t l=0; l < vec_inliers.size(); ++l)
        {
            const size_t  trackId = map_bearingIdToTrackId.at(vec_inliers[l]);
            inliers_set.insert( trackId );
        }

        openMVG::tracks::STLMAPTracks map_tracksInliers;
        size_t  counter = 0;
        for( size_t idx = 0 ; idx < map_tracks.size(); ++idx )
        {
            if(  inliers_set.find( idx ) != inliers_set.end() )
            {
                map_tracksInliers[counter] = map_tracks[idx];
                ++counter;
            }
        }

        // compute initial triangulation
        SfM_Data   tiny_scene;

        // intialize poses (which are now shared by a group of images)
        tiny_scene.poses[indexRig1] = Pose3(Mat3::Identity(), Vec3::Zero());
        tiny_scene.poses[indexRig2] = relativePose;

        // insert all views in the scene
        for (const auto & pairIterator : match_pairs )
        {
          // initialize camera indexes
          const IndexT   I = pairIterator.first;
          const IndexT   J = pairIterator.second;

          // initialize views
          const View * view_I = _sfm_data.views[I].get();
          const View * view_J = _sfm_data.views[J].get();

          // add views
          tiny_scene.views.insert(*_sfm_data.GetViews().find(pairIterator.first));
          tiny_scene.views.insert(*_sfm_data.GetViews().find(pairIterator.second));

          // add intrinsics
          tiny_scene.intrinsics.insert(*_sfm_data.GetIntrinsics().find(view_I->id_intrinsic));
          tiny_scene.intrinsics.insert(*_sfm_data.GetIntrinsics().find(view_J->id_intrinsic));
        }

        // Fill sfm_data with the computed tracks (no 3D yet)
        Landmarks & structure = tiny_scene.structure;
        IndexT idx(0);
        for (STLMAPTracks::const_iterator itTracks = map_tracksInliers.begin();
          itTracks != map_tracksInliers.end();
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

        // Compute 3D position of the landmark of the structure by triangulation of the observations
        {
          SfM_Data_Structure_Computation_Blind structure_estimator(false);
          structure_estimator.triangulate(tiny_scene);
        }

        // - refine only Structure and Rotations & translations (keep intrinsic constant)
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
        }

        {
          vec_relatives[relative_pose_pair] = std::make_pair(
            relativePose.rotation(),
            relativePose.translation());

          // TODO: Extend later to a pose graph and not a view graph
          // convert relative view index to relative pose indexes
          //const Pair relativePoseIndexes(view_I->id_pose, view_J->id_pose);
          //vec_relatives[relativePoseIndexes] = std::make_pair(
          //  relativePose.relativePose.rotation(),
          //  relativePose.relativePose.translation());
        }
      }
    }
  }
  // Log input graph to the HTML report
  if (!_sLoggingFile.empty() && !_sOutDirectory.empty())
  {
    std::set<IndexT> set_ViewIds;
      std::transform(_sfm_data.GetViews().begin(), _sfm_data.GetViews().end(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
    graph::indexedGraph putativeGraph(set_ViewIds, getPairs(_matches_provider->_pairWise_matches));
    graph::exportToGraphvizData(
      stlplus::create_filespec(_sOutDirectory, "input_largest_cc_relative_motions_graph"),
      putativeGraph.g);

    using namespace htmlDocument;
    std::ostringstream os;

    os << "<br>" << "input_largest_cc_relative_motions_graph" << "<br>"
       << "<img src=\""
       << stlplus::create_filespec(_sOutDirectory, "input_largest_cc_relative_motions_graph", "svg")
       << "\" height=\"600\">\n";

    _htmlDocStream->pushInfo(os.str());
  }
}

} // namespace sfm
} // namespace openMVG
