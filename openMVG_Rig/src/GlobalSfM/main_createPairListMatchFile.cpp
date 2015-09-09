/*
 * main_createPairListMatchFile
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

#include <cstdlib>

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/system/timer.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <unordered_map>

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Create contiguous pose pair links for view based matching:  \n"
    << "----\n"
    << " it creates a pair_list.txt file for main_ComputeMatches\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string s_SfM_Data_filename;
  std::string s_out_file;
  int i_overlapping = 5;

  cmd.add( make_option('i', s_SfM_Data_filename, "input_file") );
  cmd.add( make_option('o', s_out_file, "out_file") );
  cmd.add( make_option('p', i_overlapping, "pose_overlapping") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-o|--out_file] the output pairlist file\n"
    << "optional:\n"
    << "[-p|--pose_overlapping] Define the pose neighborhood for matching\n"
    << "\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, s_SfM_Data_filename, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< s_SfM_Data_filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (s_out_file.empty())  {
    std::cerr << "\nIt is an invalid output filename." << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(stlplus::folder_part(s_out_file)))
  {
    if (!stlplus::folder_create(stlplus::folder_part(s_out_file)))
    {
      std::cerr << "Cannot create directory for the output file." << std::endl;
      return EXIT_FAILURE;
    }
  }

  //---------------------------------------
  // List the view pose
  // establish a pose graph 'linear pose relation' with a bounded link coverage
  // conver the pose graph edges to a view graph
  //---------------------------------------
  std::unordered_multimap<IndexT, IndexT> pose_id_toViewId;
  std::set<IndexT> set_poses;

  // Get nodes of the pose graph
  for (const auto & viewIter : sfm_data.GetViews())
  {
    const View * v = viewIter.second.get();
    assert (viewIter.first == v->id_view);
    pose_id_toViewId.insert( std::make_pair(v->id_pose, v->id_view) );
    set_poses.insert(v->id_pose);
  }
  std::vector<IndexT> vec_poses(set_poses.begin(), set_poses.end());

  // Create the 'linear' pose graph pair relationship
  const Pair_Set pose_pairs = contiguousWithOverlap(set_poses.size(), i_overlapping);

  // Convert the pose graph to a view graph
  Pair_Set view_pair;
  for (const auto & pose_pair : pose_pairs)
  {
    const IndexT poseA = pose_pair.first;
    const IndexT poseB = pose_pair.second;
    // get back the view related to those poses and create the pair (exhaustively)
    const auto range_a = pose_id_toViewId.equal_range(vec_poses[poseA]);
    for (auto view_id_a = range_a.first; view_id_a != range_a.second; view_id_a++)
    {
      const auto range_b = pose_id_toViewId.equal_range(vec_poses[poseB]);
      for (auto view_id_b = range_b.first; view_id_b != range_b.second; view_id_b++)
      {
        if (view_id_a == view_id_b)
        {
          continue;
        }
        view_pair.insert(
          Pair(std::min(view_id_a->second,view_id_b->second),
               std::max(view_id_a->second,view_id_b->second)));
      }
    }
  }

  if (view_pair.empty())
  {
    std::cout << "Warning: The computed pair list is empty...!" << std::endl;
  }

  if (savePairs(s_out_file, view_pair))
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}
