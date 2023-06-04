##############################
# Preamble
##############################
cmake_minimum_required(VERSION 3.2)
list(APPEND CMAKE_MESSAGE_CONTEXT tinympc)
project(TinyMPC VERSION 0.1 LANGUAGES C CXX)

##############################
# Project wide
##############################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS NO)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS NO)

include(CMakePrintHelpers)
include(GNUInstallDirs)
# include(Functions)
# include(FetchContent)

# Add CPM Dependency Manager
# include(FindCPM)

# Handle default build type
set(TINY_DEFAULT_BUILD_TYPE "Debug")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Setting CMAKE_BUILD_TYPE to ${TINY_DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE ${TINY_DEFAULT_BUILD_TYPE} CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

# Enable testing
option(TINY_BUILD_TESTS "Build tests for TinyMPC" ON)

# Enable clang-tidy analysis
option(TINY_CLANG_TIDY "Run clang-tidy analyzer on the source code." OFF)

# Build examples
option(TINY_BUILD_EXAMPLES "Build examples for TinyMPC." ON)

# Build with -march=native
option(TINY_VECTORIZE "Compile with -march=native" ON)

# Build with extra information for profiling
option(TINY_PROFILE "Compile with -pg" OFF)

##############################
# Dependencies
##############################
include_directories(ext/Eigen)

##############################
# Main build targets
##############################
# Compile options
if (NOT WIN32)
  add_compile_options()
  # add_compile_options(-Wall -Wextra -pedantic -Wno-comment)
  # add_compile_options(-Wno-error=unknown-pragmas
  #                     -Wformat=2 -Wno-unused-parameter -Wshadow 
  #                     -Wwrite-strings -Wstrict-prototypes -Wold-style-definition 
  #                     -Wredundant-decls -Wnested-externs -Wmissing-include-dirs)
  if (TINY_VECTORIZE)
    add_compile_options(-march=native)
  endif()
  # if (TINY_PROFILE)
  #   add_compile_options(-pd)
  # endif()
endif()

if (CMAKE_C_COMPILER_ID MATCHES "GNU")
  add_compile_options(-Wlogical-op)
endif()  

# Make all includes relative to src/ directory
include_directories(${PROJECT_SOURCE_DIR}/src)

# Build source files
add_subdirectory(src/tinympc)

##############################
# Tests
##############################
if (TINY_BUILD_TESTS)
  add_subdirectory(tests)
endif()                 

##############################
# Examples
##############################
if (TINY_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()