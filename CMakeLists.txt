cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
set(PackageName LUXEDetector)
project(${PackageName})
find_package(DD4hep REQUIRED COMPONENTS DDG4)
# our sources
set(sources ./src/MLTracker_geo.cpp ./src/LUXETracker_geo.cpp)
# create our library and make the components file
add_dd4hep_plugin(${PackageName} SHARED ${sources}) # link it with DDCore, or whatever you need
target_link_libraries(${PackageName} DD4hep::DDCore)
# Create this_package.sh file, and install
dd4hep_instantiate_package(${PackageName})
