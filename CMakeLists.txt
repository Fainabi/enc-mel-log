cmake_minimum_required(VERSION 3.16)
project(ppaudio_he_frontend CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "build type" FORCE)
endif()

set(OPENFHE_SOURCE_DIR "" CACHE PATH
    "OpenFHE source tree (use with OPENFHE_BUILD_DIR for an uninstalled build)")
set(OPENFHE_BUILD_DIR "" CACHE PATH
    "OpenFHE build tree containing libOPENFHE*.so")

if(OPENFHE_SOURCE_DIR AND OPENFHE_BUILD_DIR)
  set(OpenFHE_INCLUDE
      ${OPENFHE_SOURCE_DIR}/src/pke/include
      ${OPENFHE_SOURCE_DIR}/src/core/include
      ${OPENFHE_SOURCE_DIR}/src/binfhe/include)
  set(OpenFHE_INCLUDE_PARENT ${OPENFHE_SOURCE_DIR}/src)
  set(OpenFHE_LIBDIR ${OPENFHE_BUILD_DIR}/lib)
  set(OpenFHE_SHARED_LIBRARIES OPENFHEpke OPENFHEcore OPENFHEbinfhe)
  set(BASE_OPENFHE_VERSION "source-build")
else()
  find_package(OpenFHE REQUIRED)
endif()

# OpenFHEConfig.cmake exports OpenFHE_INCLUDE as the directory that *contains*
# core/ pke/ binfhe/ (typically <prefix>/include/openfhe).  Two include roots are
# needed because the sources use both spellings of the same header:
#   src/cpp_complex_e2e.cpp : #include "math/chebyshev.h"
#   src/cpp_real_dct.cpp    : #include "openfhe/core/math/chebyshev.h"
# The first resolves against ${OpenFHE_INCLUDE}/core, the second against the
# parent of ${OpenFHE_INCLUDE}.
get_filename_component(OpenFHE_INCLUDE_PARENT "${OpenFHE_INCLUDE}" DIRECTORY)

include_directories(
  ${OpenFHE_INCLUDE}
  ${OpenFHE_INCLUDE}/core
  ${OpenFHE_INCLUDE}/pke
  ${OpenFHE_INCLUDE}/binfhe
  ${OpenFHE_INCLUDE_PARENT}
  ${OPENFHE_BUILD_DIR}/src/core
  ${OPENFHE_SOURCE_DIR}/third-party/cereal/include
)
link_directories(${OpenFHE_LIBDIR})

separate_arguments(OPENFHE_FLAGS UNIX_COMMAND "${OpenFHE_CXX_FLAGS}")

set(PPAUDIO_TARGETS
  cpp_complex_e2e
  cpp_complex_ops
  cpp_complex_smoke
  cpp_dct_small
  cpp_inv2_placement
  cpp_linear_smoke
  cpp_mv_unit
  cpp_real_dct
)

foreach(tgt IN LISTS PPAUDIO_TARGETS)
  add_executable(${tgt} src/${tgt}.cpp)
  target_compile_options(${tgt} PRIVATE ${OPENFHE_FLAGS})
  target_link_libraries(${tgt} PRIVATE ${OpenFHE_SHARED_LIBRARIES})
endforeach()

message(STATUS "OpenFHE include : ${OpenFHE_INCLUDE}")
message(STATUS "OpenFHE libdir  : ${OpenFHE_LIBDIR}")
message(STATUS "OpenFHE version : ${BASE_OPENFHE_VERSION}")
message(STATUS "targets         : ${PPAUDIO_TARGETS}")
