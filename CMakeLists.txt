# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(Software-Development-for-Algorithmic-Problems)


# CGAL and its components
find_package(CGAL QUIET COMPONENTS Qt5)

if (NOT CGAL_FOUND)

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package(Boost REQUIRED)

if (NOT Boost_FOUND)

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories(BEFORE includes)

# include for local package


# Creating entries for all C++ files with "main" routine
# ##########################################################
create_single_source_cgal_program("main.cpp")

add_executable(main src/main.cpp src/parsing.cpp src/dtriangulation.cpp src/edgeflip.cpp)

if (CGAL_Qt5_FOUND)
    add_definitions(-DCGAL_USE_BASIC_VIEWER)
    target_link_libraries(main PUBLIC CGAL::CGAL_Qt5)
endif()