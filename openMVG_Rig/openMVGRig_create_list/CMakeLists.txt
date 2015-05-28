 #
 #  listomvg
 #
 # Copyright (c) 2013-2014 FOXEL SA - http://foxel.ch
 # Please read <http://foxel.ch/license> for more information.
 #
 #
 # Author(s):
 #
 #      Stéphane Flotron <s.flotron@foxel.ch>
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


cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# ==============================================================================
# Build Executable
# ==============================================================================
add_executable(listomvg
                    tools.cpp
                    types.cpp
                    list_utils.cpp
                    list.cpp
                    )
target_link_libraries(listomvg
  ${LISTOMVG_LIBRARY_LIST}
)

add_dependencies(listomvg fastcal stlplus)
set_target_properties( listomvg PROPERTIES RUNTIME_OUTPUT_DIRECTORY .. )

install(TARGETS listomvg DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
