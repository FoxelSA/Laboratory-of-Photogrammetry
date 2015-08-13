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
#include "third_party/htmlDoc/htmlDoc.hpp"

#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/multiview/essential.hpp"

#include "third_party/progress/progress.hpp"

#include <opengv/types.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>

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

  for (const auto & relativePoseIterator : rigWiseMatches)
  {
    const Pair relative_pose_pair = relativePoseIterator.first;
    const Pair_Set & match_pairs = relativePoseIterator.second;

    // Compute the relative pose...
    // - if only one pair of matches: relative pose between two pinhole images
    // - if many pairs: relative pose between rigid rigs
  }
}

} // namespace sfm
} // namespace openMVG
