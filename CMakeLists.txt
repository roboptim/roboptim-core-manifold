# Copyright 2015, Benjamin Chrétien, CNRS-LIRMM.
#
# This file is part of roboptim-core.
# roboptim-core is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# roboptim-core is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with roboptim-core.  If not, see <http://www.gnu.org/licenses/>.

# Requires at least CMake 2.8 to configure the package.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/lapack.cmake)
INCLUDE(cmake/cpack.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME roboptim-core-manifold)
SET(PROJECT_DESCRIPTION "RobOptim Core - manifold support")
SET(PROJECT_URL "http://github.com/roboptim/roboptim-core-manifold")
SET(PROJECT_DEBUG_POSTFIX "_d")

SET(HEADERS
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/descriptive-wrapper.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/descriptive-wrapper.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/function-on-manifold.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/function-on-manifold.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-desc.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-desc.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-map.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/problem-on-manifold.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/problem-on-manifold.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-problem-factory.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-problem-factory.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-merger.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/wrapper-on-manifold.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/wrapper-on-manifold.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/sum-on-manifold.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/sum-on-manifold.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/util.hh
  )

SETUP_PROJECT()

# If compiler support symbol visibility, enable it.
INCLUDE(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG(-fvisibility=hidden HAS_VISIBILITY)
IF(HAS_VISIBILITY)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
ENDIF()
CHECK_C_COMPILER_FLAG(-ferror-limit=1 HAS_ERROR_LIMIT)
IF(HAS_ERROR_LIMIT)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ferror-limit=1")
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O1 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")


# Search for dependencies.
SET(BOOST_COMPONENTS
  date_time filesystem system thread program_options unit_test_framework)
SEARCH_FOR_BOOST()

ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.2.0")
ADD_REQUIRED_DEPENDENCY("roboptim-core >= 3.1")
ADD_REQUIRED_DEPENDENCY("manifolds" "manifolds_d")

HEADER_INSTALL("${HEADERS}")

OPTION(DISABLE_TESTS "Disable test programs" OFF)

ADD_SUBDIRECTORY(src)

IF(NOT DISABLE_TESTS)
  ADD_SUBDIRECTORY(tests)
ELSE()
  MESSAGE(WARNING
    "Tests should only be disabled for speficic cases. Do it at your own risk.")
ENDIF()

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()
