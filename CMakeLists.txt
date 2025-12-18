cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(luxegeo)

# install(DIRECTORY compact DESTINATION share/luxegeo/)
# SET( CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" )
SET( CMAKE_INSTALL_LIBDIR "lib" )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

SET( ${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}" )

find_package(DD4hep REQUIRED COMPONENTS DDG4)
# create our library and make the components file
# our sources
set(sources
  ./src/SimpleTracker_geo.cpp
  ./src/LUXETracker_geo.cpp
  ./src/LUXEEcal_geo.cpp
  )

add_dd4hep_plugin(${PROJECT_NAME} SHARED ${sources}) # link it with DDCore, or whatever you need
target_link_libraries(${PROJECT_NAME} DD4hep::DDCore)

# Create this_package.sh file, and install
dd4hep_instantiate_package(${PROJECT_NAME})

# install(DIRECTORY compact DESTINATION share/luxegeo/)

include(CTest)
add_subdirectory(tests)
