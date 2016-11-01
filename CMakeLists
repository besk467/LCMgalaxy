cmake_minimum_required (VERSION 2.6)

project (Galaxy)

# global property to hold the names of katydid library targets
set_property (GLOBAL PROPERTY ${PROJECT_NAME}_LIBRARIES)

include_directories(/usr/local/installed_packages/root/lib)

# deal with the rpath settings (from http://www.cmake.org/Wiki/CMake_RPATH_handling)
# use (i.e. don't skip) the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_BUILD_DIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

# add to the RPATH to be used when installing, but only if it's not a system directory
#list (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${LIB_INSTALL_DIR}" isSystemDir)
#if ("${isSystemDir}" STREQUAL "-1")
#   list (APPEND CMAKE_INSTALL_RPATH "${LIB_INSTALL_DIR}")
#endif ("${isSystemDir}" STREQUAL "-1")


set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}")

include_directories (${PROJECT_SOURCE_DIR})

find_package (ROOT 5.28)
if (NOT ROOT_FOUND)
    message (FATAL "ROOT not found")
endif (NOT ROOT_FOUND)
include_directories (${ROOT_INCLUDES})

#libMinuit doesn't get included in ROOT_LIBRARIES for some reason
set( ROOT_LIBRARIES ${ROOT_LIBRARIES} ${ROOT_LIBRARY_DIR}/libMinuit.so)

set( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wl,-rpath,${ROOT_LIBRARY_DIR} )

set (HEADER_FILES
    Galaxy.hh
    GalaxyDB.hh
    CalcH.hh
    Constants.hh
)

set (SOURCE_FILES
    Galaxy.cc
    GalaxyDB.cc
    CalcH.cc
    KDFit3e_prog.cc
)

add_executable (KDFit3e ${SOURCE_FILES})
target_link_libraries (KDFit3e ${ROOT_LIBRARIES})
