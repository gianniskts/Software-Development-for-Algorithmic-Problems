# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(Software-Development-for-Algorithmic-Problems)

set(CMAKE_SKIP_RPATH TRUE)

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
include_directories(
  BEFORE 
  ${CMAKE_CURRENT_SOURCE_DIR}/includes
)

set(SOURCES
  main.cpp
  src/Parsing.cpp
  src/PolygonManipulation.cpp
  src/ActionFunctions.cpp
  src/MyTriangulation.cpp
)

add_executable(main ${SOURCES})

if (CGAL_Qt5_FOUND)
    target_compile_definitions(main PRIVATE CGAL_USE_BASIC_VIEWER)
    target_link_libraries(main PUBLIC CGAL::CGAL_Qt5)
endif()

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)