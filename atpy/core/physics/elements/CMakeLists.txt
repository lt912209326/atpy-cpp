aux_source_directory(. ELEMENTS_SRC_LIST)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../utils
    )

IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++14 -fPIC  -fopenmp -ffast-math")
    SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../Release)
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "Windows" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /std:c++14  /openmp")
    SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../..)
ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")

add_library(Elements  STATIC ${ELEMENTS_SRC_LIST} )
