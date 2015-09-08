/*
 * main_openMVG2DAV
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
#include "openMVG/system/timer.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"


int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Convert rig poses to txt pose file:   \n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "Cannot create the output directory" << std::endl;
      return EXIT_FAILURE;
    }

  //---------------------------------------
  // List the view and log the global pose
  //---------------------------------------
  std::set<IndexT> handled_poses;

  C_Progress_display progress_bar( sfm_data.GetViews().size() );
  for (const auto & viewIter : sfm_data.GetViews())
  {
    ++progress_bar;
    const View * v = viewIter.second.get();
    if (!sfm_data.IsPoseAndIntrinsicDefined(v))
    {
      continue;
    }
    if (handled_poses.insert(v->id_pose).second)
    {
      // extract timestamp
      const std::string viewName = v->s_Img_path;
      std::size_t delimiter = viewName.find_first_of("-");
      if (delimiter == std::string::npos)
      {
        continue;
      }
      // export the pose
      std::ofstream ofs(stlplus::create_filespec(sOutDir, viewName.substr (0,delimiter), "txt").c_str());
      const Mat3 R = sfm_data.poses.at(v->id_pose).rotation();
      const Vec C = sfm_data.poses.at(v->id_pose).center();
      ofs
        << R(0,0) << " " << R(0,1) << " " << R(0,2) << "\n"
        << R(1,0) << " " << R(1,1) << " " << R(1,2) << "\n"
        << R(2,0) << " " << R(2,1) << " " << R(2,2) << "\n"
        << C(0)   << " " << C(1)   << " " << C(2);
      ofs.close();
    }
  }
  return EXIT_SUCCESS;
}
