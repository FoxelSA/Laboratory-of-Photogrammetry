 #
 #  OPENMVG_RIG
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


cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(openMVG_Rig C CXX)
add_definitions(-std=c++11)

# ==============================================================================
#  C / C++ flags
# ==============================================================================
set(CMAKE_CXX_COMPILER c++)
set(CMAKE_BUILD_TYPE Debug)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(warnings "-Wall -Wextra -pedantic")
    set(cxx_flags "-std=c++11")
endif()

if (NOT CONFIGURED_ONCE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings} ${cxx_flags}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif()

set(OPENMVG_RIG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )

# ==============================================================================
# OpenMP detection
# ==============================================================================
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  if (NOT MSVC)
    LIST(APPEND OPENMVG_RIG_LIBRARY_DEPENDENCIES gomp)
  endif (NOT MSVC)
endif (OPENMP_FOUND)


# ==============================================================================
# OPENMVG_RIG check for dependencies
# ==============================================================================

# ------------------------------------------------------------------------------
# pkgconfig / opencv
# ------------------------------------------------------------------------------
find_package(PkgConfig)
find_package(OpenMVG REQUIRED)

# ------------------------------------------------------------------------------
# libfastcal
# ------------------------------------------------------------------------------
add_subdirectory(lib)
set(LIBFASTCAL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/libfastcal/src/)
set(LIBFASTCAL_LIBS ${CMAKE_CURRENT_BINARY_DIR}/lib/libfastcal/liblibfastcal.a)

# ==============================================================================
# Include directories
# ==============================================================================
include_directories(
  ${OPENMVG_RIG_SOURCE_DIR}/openMVGRig_create_list
  ${OPENMVG_INCLUDE_DIRS}
  ${LIBFASTCAL_INCLUDE_DIR}
)

# ==============================================================================
# Library list
# ==============================================================================
set(OPENMVG_RIG_LIBRARY_LIST
  ${LIBFASTCAL_LIBS}
  ${OPENMVG_LIBRARIES}
)

# ==============================================================================
# Go to subdirectory
# ==============================================================================
add_subdirectory(openMVGRig_create_list)
