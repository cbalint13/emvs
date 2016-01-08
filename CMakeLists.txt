#/*
# *  Copyright (c) 2014  Balint Cristian (cristian.balint@gmail.com)
# *
# * Redistribution and use in source and binary forms, with or without
# * modification, are permitted provided that the following conditions are met:
# *
# * * Redistributions of source code must retain the above copyright notice, this
# *   list of conditions and the following disclaimer.
# *
# * * Redistributions in binary form must reproduce the above copyright notice,
# *   this list of conditions and the following disclaimer in the documentation
# *   and/or other materials provided with the distribution.
# *
# * * Neither the name of copyright holder nor the names of its
# *   contributors may be used to endorse or promote products derived from
# *   this software without specific prior written permission.
# *
# * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# *
# */
#
#/* CMakeLists.txt */
#/* Project builder */

PROJECT(emvs)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(COMMAND CMAKE_POLICY)
     CMAKE_POLICY(SET CMP0016 NEW)
ENDIF(COMMAND CMAKE_POLICY)

SET(CMAKE_SKIP_BUILD_RPATH TRUE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

SET(CMAKE_MODULE_PATH "cmake/;${CMAKE_MODULE_PATH}")

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)


# Set a default build type if none was specified
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  MESSAGE(STATUS "Setting build type to 'Release' as none was specified.")
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
ENDIF()

SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3 -msse2 -mssse3 -msse4.2 -ffast-math -funroll-loops")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++11 -g -O0")

#####
# find OpenMP library
FIND_PACKAGE(OpenMP QUIET)
IF(OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  MESSAGE(STATUS "OpenMP found.")
ENDIF()

#####
# find OpenCV
FIND_PACKAGE(OpenCV REQUIRED core highgui imgproc imgcodecs xfeatures2d)
IF(OpenCV_FOUND)
  INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
  MESSAGE(STATUS "OpenCV version: " ${OpenCV_VERSION})
  STRING(COMPARE LESS "${OpenCV_VERSION}" 3.0.0 VersionIncompatible)
  IF(VersionIncompatible)
    MESSAGE(FATAL_ERROR "OpenCV minimum 3.0.0 required.")
  ENDIF()
ENDIF()

# local directory
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include/)

ADD_SUBDIRECTORY(src)
