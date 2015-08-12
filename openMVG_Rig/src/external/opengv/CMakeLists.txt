cmake_minimum_required(VERSION 2.4.6)
#project(opengv)

set (OPENGV_VERSION_MAJOR 1)
set (OPENGV_VERSION_MINOR 0)

# Set the build type.  Options are:
#
#  None (CMAKE_C_FLAGS or CMAKE_CXX_FLAGS used)
#  Debug (CMAKE_C_FLAGS_DEBUG or CMAKE_CXX_FLAGS_DEBUG)
#  Release (CMAKE_C_FLAGS_RELEASE or CMAKE_CXX_FLAGS_RELEASE)
#  RelWithDebInfo (CMAKE_C_FLAGS_RELWITHDEBINFO or CMAKE_CXX_FLAGS_RELWITHDEBINFO
#  MinSizeRel (CMAKE_C_FLAGS_MINSIZEREL or CMAKE_CXX_FLAGS_MINSIZEREL)

set(CMAKE_BUILD_TYPE Release)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# get eigen (and boost under windows)
set(OPENGV_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
set(CMAKE_MODULE_PATH_OLD ${CMAKE_MODULE_PATH})

IF(WIN32)
  set(ADDITIONAL_INCLUDE_DIRS ${OPENGV_SOURCE_DIR}/third_party
                              ${OPENGV_SOURCE_DIR}/third_party/eigen3
                              ${OPENGV_SOURCE_DIR}/third_party/eigen3/unsupported )
ELSE()
  unset( CMAKE_MODULE_PATH )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${OPENGV_SOURCE_DIR}/modules/")
  find_package(Eigen REQUIRED)
  set(ADDITIONAL_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIR}/unsupported )
ENDIF()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_OLD})

include_directories(${OPENMVG_INCLUDE_DIRS} ${OPENGV_SOURCE_DIR}/include)

set( OPENGV_SOURCE_FILES
  src/absolute_pose/modules/main.cpp
  src/absolute_pose/modules/gp3p/code.cpp
  src/absolute_pose/modules/gp3p/init.cpp
  src/absolute_pose/modules/gp3p/reductors.cpp
  src/absolute_pose/modules/gp3p/spolynomials.cpp
  src/absolute_pose/modules/Epnp.cpp
  src/absolute_pose/modules/gpnp1/code.cpp
  src/absolute_pose/modules/gpnp1/init.cpp
  src/absolute_pose/modules/gpnp1/reductors.cpp
  src/absolute_pose/modules/gpnp1/spolynomials.cpp
  src/absolute_pose/modules/gpnp2/code.cpp
  src/absolute_pose/modules/gpnp2/init.cpp
  src/absolute_pose/modules/gpnp2/reductors.cpp
  src/absolute_pose/modules/gpnp2/spolynomials.cpp
  src/absolute_pose/modules/gpnp3/code.cpp
  src/absolute_pose/modules/gpnp3/init.cpp
  src/absolute_pose/modules/gpnp3/reductors.cpp
  src/absolute_pose/modules/gpnp3/spolynomials.cpp
  src/absolute_pose/modules/gpnp4/code.cpp
  src/absolute_pose/modules/gpnp4/init.cpp
  src/absolute_pose/modules/gpnp4/reductors.cpp
  src/absolute_pose/modules/gpnp4/spolynomials.cpp
  src/absolute_pose/modules/gpnp5/code.cpp
  src/absolute_pose/modules/gpnp5/init.cpp
  src/absolute_pose/modules/gpnp5/reductors.cpp
  src/absolute_pose/modules/gpnp5/spolynomials.cpp
  src/absolute_pose/modules/upnp2.cpp
  src/absolute_pose/modules/upnp4.cpp
  src/relative_pose/modules/main.cpp
  src/relative_pose/modules/fivept_nister/modules.cpp
  src/relative_pose/modules/fivept_stewenius/modules.cpp
  src/relative_pose/modules/fivept_kneip/code.cpp
  src/relative_pose/modules/fivept_kneip/init.cpp
  src/relative_pose/modules/fivept_kneip/reductors.cpp
  src/relative_pose/modules/fivept_kneip/spolynomials.cpp
  src/relative_pose/modules/sixpt/modules2.cpp
  src/relative_pose/modules/eigensolver/modules.cpp
  src/relative_pose/modules/ge/modules.cpp
  src/math/cayley.cpp
  src/math/quaternion.cpp
  src/math/arun.cpp
  src/math/Sturm.cpp
  src/math/roots.cpp
  src/math/gauss_jordan.cpp
  src/absolute_pose/methods.cpp
  src/absolute_pose/CentralAbsoluteAdapter.cpp
  src/absolute_pose/NoncentralAbsoluteAdapter.cpp
  src/relative_pose/methods.cpp
  src/relative_pose/CentralRelativeAdapter.cpp
  src/relative_pose/CentralRelativeWeightingAdapter.cpp
  src/relative_pose/NoncentralRelativeAdapter.cpp
  src/relative_pose/CentralRelativeMultiAdapter.cpp
  src/relative_pose/NoncentralRelativeMultiAdapter.cpp
  src/triangulation/methods.cpp
  src/point_cloud/methods.cpp
  src/point_cloud/PointCloudAdapter.cpp
  src/sac_problems/absolute_pose/AbsolutePoseSacProblem.cpp
  src/sac_problems/relative_pose/CentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/RotationOnlySacProblem.cpp
  src/sac_problems/relative_pose/TranslationOnlySacProblem.cpp
  src/sac_problems/relative_pose/EigensolverSacProblem.cpp
  src/sac_problems/relative_pose/MultiCentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.cpp
  src/sac_problems/point_cloud/PointCloudSacProblem.cpp
  src/absolute_pose/MACentralAbsolute.cpp
  src/absolute_pose/MANoncentralAbsolute.cpp
  src/relative_pose/MACentralRelative.cpp
  src/relative_pose/MANoncentralRelative.cpp
  src/relative_pose/MANoncentralRelativeMulti.cpp
  src/point_cloud/MAPointCloud.cpp )

set( OPENGV_HEADER_FILES
  include/opengv/types.hpp
  include/opengv/OptimizationFunctor.hpp
  include/opengv/absolute_pose/methods.hpp
  include/opengv/relative_pose/methods.hpp
  include/opengv/triangulation/methods.hpp
  include/opengv/point_cloud/methods.hpp
  include/opengv/math/cayley.hpp
  include/opengv/math/quaternion.hpp
  include/opengv/math/arun.hpp
  include/opengv/math/Sturm.hpp
  include/opengv/math/roots.hpp
  include/opengv/math/gauss_jordan.hpp
  include/opengv/absolute_pose/AbsoluteAdapterBase.hpp
  include/opengv/absolute_pose/CentralAbsoluteAdapter.hpp
  include/opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp
  include/opengv/absolute_pose/AbsoluteMultiAdapterBase.hpp
  include/opengv/relative_pose/RelativeAdapterBase.hpp
  include/opengv/relative_pose/RelativeMultiAdapterBase.hpp
  include/opengv/relative_pose/CentralRelativeAdapter.hpp
  include/opengv/relative_pose/CentralRelativeWeightingAdapter.hpp
  include/opengv/relative_pose/NoncentralRelativeAdapter.hpp
  include/opengv/relative_pose/CentralRelativeMultiAdapter.hpp
  include/opengv/relative_pose/NoncentralRelativeMultiAdapter.hpp
  include/opengv/point_cloud/PointCloudAdapterBase.hpp
  include/opengv/point_cloud/PointCloudAdapter.hpp
  include/opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/MultiCentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/EigensolverSacProblem.hpp
  include/opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp
  include/opengv/sac_problems/relative_pose/TranslationOnlySacProblem.hpp
  include/opengv/sac_problems/point_cloud/PointCloudSacProblem.hpp
  include/opengv/absolute_pose/MACentralAbsolute.hpp
  include/opengv/absolute_pose/MANoncentralAbsolute.hpp
  include/opengv/relative_pose/MACentralRelative.hpp
  include/opengv/relative_pose/MANoncentralRelative.hpp
  include/opengv/relative_pose/MANoncentralRelativeMulti.hpp
  include/opengv/point_cloud/MAPointCloud.hpp )

IF(WIN32)
  add_library( opengv ${OPENGV_SOURCE_FILES} ${OPENGV_HEADER_FILES} )
  add_library( random_generators test/random_generators.cpp test/random_generators.hpp test/experiment_helpers.cpp test/experiment_helpers.hpp test/time_measurement.cpp test/time_measurement.hpp )
ELSE()
  add_library( opengv SHARED ${OPENGV_SOURCE_FILES} ${OPENGV_HEADER_FILES} )
  add_library( random_generators SHARED test/random_generators.cpp test/random_generators.hpp test/experiment_helpers.cpp test/experiment_helpers.hpp test/time_measurement.cpp test/time_measurement.hpp )
  target_link_libraries(random_generators opengv)
ENDIF()
