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

SET(HEADERS
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/descriptive-wrapper.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/descriptive-wrapper.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/dispatcher.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/function-on-manifold.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/function-on-manifold.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-desc.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-desc.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/manifold-map.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/problem-on-manifold.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/problem-on-manifold.hxx
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/problem-factory.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/manifold-map/decorator/problem-factory.hxx
  )

SETUP_PROJECT()

# If compiler support symbol visibility, enable it.
INCLUDE(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG(-fvisibility=hidden HAS_VISIBILITY)
IF(HAS_VISIBILITY)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -std=c++0x")
ENDIF()

# TODO: Add main library to pkg-config file.
#PKG_CONFIG_APPEND_LIBS(roboptim-core-manifold)

# Search for dependencies.
SET(BOOST_COMPONENTS
  date_time filesystem system thread program_options unit_test_framework)
SEARCH_FOR_BOOST()

ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.2.0")
ADD_REQUIRED_DEPENDENCY("liblog4cxx >= 0.10.0")
ADD_REQUIRED_DEPENDENCY("roboptim-core >= 2.0")

ADD_REQUIRED_DEPENDENCY("manifolds")

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
