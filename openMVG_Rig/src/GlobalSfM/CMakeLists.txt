 #
 #  OPENMVG_RIG
 #
 # Copyright (c) 2013-2014 FOXEL SA - http://foxel.ch
 # Please read <http://foxel.ch/license> for more information.
 #
 #
 # Author(s):
 #
 #      Pierre Moulon <p.moulon@foxel.ch>
 #
 #
 # This file is part of the FOXEL project <http://foxel.ch>.
 #
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Affero General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU Affero General Public License for more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with this program.  if not, see <http://www.gnu.org/licenses/>.
 #
 #
 # Additional Terms:
 #
 #      You are required to preserve legal notices and author attributions in
 #      that material or in the Appropriate Legal Notices displayed by works
 #      containing it.
 #
 #      You are required to attribute the work as explained in the "Usage and
 #      Attribution" section of <http://foxel.ch/license>.
 #/

# ==============================================================================
# Build Executable
# ==============================================================================
add_executable(openMVGRig_GlobalSfM
  sfm_engine_rigid_rig.hpp
  sfm_engine_rigid_rig.cpp
  main_GlobalSfM.cpp)

add_dependencies(openMVGRig_GlobalSfM opengv )
target_link_libraries(openMVGRig_GlobalSfM ${OPENMVG_LIBRARIES})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/software)

install(TARGETS openMVGRig_GlobalSfM DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
