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

#include "openMVG/sfm/pipelines/sfm_engine.hpp"

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

namespace openMVG{
namespace sfm{

/// A SfM Pipeline Reconstruction Engine that is adapted for:
///  - calibrated camera rig (known intrinsics & subposes).
/// Method:
/// - global fusion of relative motions.
class ReconstructionEngine_RelativeMotions_RigidRig : public ReconstructionEngine
{
public:

  ReconstructionEngine_RelativeMotions_RigidRig(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~ReconstructionEngine_RelativeMotions_RigidRig();

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);

  void SetRotationAveragingMethod(ERotationAveragingMethod eRotationAveragingMethod);
  void SetTranslationAveragingMethod(ETranslationAveragingMethod _eTranslationAveragingMethod);

  virtual bool Process();

protected:
  /// Compute from relative rotations the global rotations of the camera poses
  bool Compute_Global_Rotations
  (
    const RelativeInfo_Map & vec_relatives,
    Hash_Map<IndexT, Mat3> & map_globalR
  );

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations
  );

  /// Compute the initial structure of the scene
  bool Compute_Initial_Structure();

  // Adjust the scene (& remove outliers)
  bool Adjust();

private:
  /// Compute relative rotations
  void Compute_Relative_Rotations(RelativeInfo_Map & vec_relatives);

  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> _htmlDocStream;
  std::string _sLoggingFile;

  // Parameter
  ERotationAveragingMethod _eRotationAveragingMethod;
  ETranslationAveragingMethod _eTranslationAveragingMethod;

  //-- Data provider
  Features_Provider  * _features_provider;
  Matches_Provider  * _matches_provider;

  std::shared_ptr<Features_Provider> _normalized_features_provider;
};

} // namespace sfm
} // namespace openMVG
