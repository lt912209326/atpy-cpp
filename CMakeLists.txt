cmake_minimum_required(VERSION 3.10)
project(atpy_cpp)



find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  message(STATUS "Set up ccache ...")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif()

aux_source_directory(. ATPY_CPP_SRC_LIST)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/atpy
    ${CMAKE_CURRENT_SOURCE_DIR}/atpy/core
    ${CMAKE_CURRENT_SOURCE_DIR}/atpy/core/physics
    ${CMAKE_CURRENT_SOURCE_DIR}/atpy/core/physics/beamline
    ${CMAKE_CURRENT_SOURCE_DIR}/atpy/core/physics/elements
    ${CMAKE_CURRENT_SOURCE_DIR}/atpy/core/physics/utils
    )
# set(CMAKE_AUTOMOC ON)


set(CMAKE_CXX_STANDARD 14)
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++14 -fPIC -fopenmp -ffast-math")
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "Windows" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /std:c++14  /openmp")
ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/atpy/core/physics/elements)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/atpy/core/physics/beamline)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/atpy/core)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /NODEFAULTLIB:libcmt.lib /std:c++17  /openmp")

# target_link_libraries(
#     test_acc 
#     # ${CMAKE_CURRENT_LIST_DIR}/src/physics/lib/beamline/
#     BeamLine
#     Elements
# )
