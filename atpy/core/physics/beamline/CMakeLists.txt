aux_source_directory(. BEAMLINE_SRC_LIST)
aux_source_directory(../utils UTILS_SRC_LIST)
# aux_source_directory(../elements ELEMENTS_SRC_LIST)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../Eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/../elements
    ${CMAKE_CURRENT_SOURCE_DIR}/../utils
    )


IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++14 -fPIC  -fopenmp -ffast-math")
    SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../Release)
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "Windows" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /std:c++14  /openmp")
    SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../..)
ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..)
# add_library(BeamLine  ${BEAMLINE_SRC_LIST} ${UTILS_SRC_LIST} ${ELEMENTS_SRC_LIST} )


add_library(BeamLine STATIC ${BEAMLINE_SRC_LIST} ${UTILS_SRC_LIST} )
