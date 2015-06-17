/*
 * main.cpp
 *
 * Copyright (c) 2014-2015 FOXEL SA - http://foxel.ch
 * Please read <http://foxel.ch/license> for more information.
 *
 *
 * Author(s):
 *
 *      Pierre Moulon <p.moulon@foxel.ch>
 *
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

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <map>

using namespace openMVG;
using namespace openMVG::features;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace svg;

// Convert HUE color to RGB
inline float hue2rgb(float p, float q, float t){
  if(t < 0) t += 1;
  if(t > 1) t -= 1;
  if(t < 1.f/6.f) return p + (q - p) * 6.f * t;
  if(t < 1.f/2.f) return q;
  if(t < 2.f/3.f) return p + (q - p) * (2.f/3.f - t) * 6.f;
  return p;
}

//
// Converts an HSL color value to RGB. Conversion formula
// adapted from http://en.wikipedia.org/wiki/HSL_color_space.
// Assumes h, s, and l are contained in the set [0, 1] and
// returns r, g, and b in the set [0, 255].
void hslToRgb(
  float h, float s, float l,
  unsigned char & r, unsigned char & g, unsigned char & b)
{
  if(s == 0){
    r = g = b = static_cast<unsigned char>(l * 255.f); // achromatic
  }else{
    const float q = l < 0.5f ? l * (1 + s) : l + s - l * s;
    const float p = 2.f * l - q;
    r = static_cast<unsigned char>(hue2rgb(p, q, h + 1.f/3.f) * 255.f);
    g = static_cast<unsigned char>(hue2rgb(p, q, h) * 255.f);
    b = static_cast<unsigned char>(hue2rgb(p, q, h - 1.f/3.f) * 255.f);
  }
}

int main(int argc, char ** argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sMatchFile;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('d', sMatchesDir, "matchdir") );
  cmd.add( make_option('m', sMatchFile, "matchfile") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Export rig-wise matches.\n"
        << "Usage: " << argv[0] << "\n"
        << "[-i|--input_file file] path to a SfM_Data scene\n"
        << "[-d|--matchdir path] path to match directory\n"
        << "[-m|--sMatchFile filename] match file (matches.x.txt)\n"
        << "[-o|--outdir path] path where the matches will be exported\n"
        << std::endl
        << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }


  //---------------------------------------
  // Read SfM Scene (image view names)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Read the features
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }

  // Read the matches
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if (!matches_provider->load(sfm_data, sMatchFile)) {
    std::cerr << std::endl
      << "Invalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  // ------------
  // For each pose pair, export the matches
  // ------------

  /// pairwise view relation shared between poseIds
  typedef std::map< Pair, Pair_Set > RigWiseMatches;

  // List shared correspondences between pose_pair
  RigWiseMatches rigWiseMatches;
  for (PairWiseMatches::const_iterator iterMatches = matches_provider->_pairWise_matches.begin();
    iterMatches != matches_provider->_pairWise_matches.end(); ++iterMatches)
  {
    const Pair pair = iterMatches->first;
    const View * v1 = sfm_data.GetViews().at(pair.first).get();
    const View * v2 = sfm_data.GetViews().at(pair.second).get();
    rigWiseMatches[Pair(v1->id_pose, v2->id_pose)].insert(pair);
  }

  //--
  // Display the Rigwise matches as SVG files
  //--

  stlplus::folder_create(sOutDir);
  std::cout << "\n Export rig-wise matches" << std::endl;
  C_Progress_display my_progress_bar( rigWiseMatches.size() );

  const bool b_yOffset = true;
  for (const auto & keyVal : rigWiseMatches)
  {
    ++my_progress_bar;

    const Pair pose_pair = keyVal.first;
    size_t x_offset = 0;
    svgDrawer svgStream;
    // Loop along the pairwise View relations contained in the rig
    for (const auto & value : keyVal.second)
    {
      const Pair view_pair = value;
      const View * view_I = sfm_data.GetViews().at(view_pair.first).get();
      const View * view_J = sfm_data.GetViews().at(view_pair.second).get();
      const std::string sView_I = stlplus::create_filespec(sfm_data.s_root_path, view_I->s_Img_path);
      const std::string sView_J = stlplus::create_filespec(sfm_data.s_root_path, view_J->s_Img_path);

      const IndMatches & matches = matches_provider->_pairWise_matches.at(view_pair);
      if (matches.empty())
      {
        continue;
      }

      // Draw images
      svgStream.drawImage(sView_I, view_I->ui_width, view_I->ui_height, x_offset);
      svgStream.drawImage(sView_J, view_J->ui_width, view_J->ui_height, x_offset, view_I->ui_height);

      //--
      // Draw matches
      //--

      const PointFeatures & vec_feat_I = feats_provider->getFeatures(view_I->id_view);
      const PointFeatures & vec_feat_J = feats_provider->getFeatures(view_J->id_view);

      const float yOffset = b_yOffset ? view_I->ui_height : 0;

      //-- Draw link between features :
      for (size_t i=0; i< matches.size(); ++i)  {
        const PointFeature & imaA = vec_feat_I[matches[i]._i];
        const PointFeature & imaB = vec_feat_J[matches[i]._j];
        unsigned char r, g, b;
        hslToRgb( (rand()%360) / 360., 1.0, .5, r, g, b);
        std::ostringstream osCol;
        osCol << "rgb(" << (int)r <<',' << (int)g << ',' << (int)b <<")";
        svgStream.drawLine(imaA.x() + x_offset, imaA.y(),
          imaB.x() + x_offset, imaB.y() + yOffset, svgStyle().stroke(osCol.str(), 2.0));
      }

      //-- Draw features (in two loop, in order to have the features upper the link, svg layer order):
      for (size_t i=0; i< matches.size(); ++i)  {
        const PointFeature & imaA = vec_feat_I[matches[i]._i];
        const PointFeature & imaB = vec_feat_J[matches[i]._j];
        svgStream.drawCircle(imaA.x() + x_offset, imaA.y(), 3.0,
          svgStyle().stroke("yellow", 2.0));
        svgStream.drawCircle(imaB.x() + x_offset, imaB.y() + yOffset, 3.0,
          svgStyle().stroke("yellow", 2.0));
      }

      x_offset += std::max(view_I->ui_width, view_J->ui_width);
    }
    std::ostringstream os;
    os << stlplus::folder_append_separator(sOutDir)
      << pose_pair.first << "_" << pose_pair.second
      << "_" << keyVal.second.size() << "_.svg";
    std::ofstream svgFile( os.str().c_str() );
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
  }
  return EXIT_SUCCESS;
}
