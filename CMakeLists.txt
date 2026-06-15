cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(lxelgeo)

# install(DIRECTORY compact DESTINATION share/luxegeo/)
# SET( CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" )
SET( CMAKE_INSTALL_LIBDIR "lib" )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

SET( ${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}" )

find_package(DD4hep REQUIRED COMPONENTS DDG4)

include_directories(${PROJECT_NAME} include)
# create our library and make the components file
# our sources
set(sources
  ./src/LxXTD8.cpp
  ./src/LxAux.cpp
  ./src/LxInteractionChamber.cpp
  ./src/FlashMagnetAssembly.cpp
  ./src/LxMagnet.cpp
  ./src/LxBeamPipes.cpp
  ./src/LxBeamPipes_geo.cpp
  ./src/LxDetChamberAssembly.cpp
  ./src/LxHICSBeamDump.cpp
  )

add_dd4hep_plugin(${PROJECT_NAME} SHARED ${sources}) # link it with DDCore, or whatever you need
target_link_libraries(${PROJECT_NAME} DD4hep::DDCore)

# Create this_package.sh file, and install
dd4hep_instantiate_package(${PROJECT_NAME})

install(DIRECTORY compact DESTINATION share/lxelgeo/)

#include(CTest)
#add_subdirectory(tests)

target_compile_options(${PROJECT_NAME} PUBLIC -ggdb -g)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ggdb -g")

