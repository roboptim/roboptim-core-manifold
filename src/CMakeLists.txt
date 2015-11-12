# Copyright 2015, Benjamin Chrétien, CNRS-LIRMM
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

# Note: header-only, so no library to compile here.

ADD_LIBRARY(roboptim-core-manifold SHARED
  ${HEADERS}
  manifold-merger.cc)

IF(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILDING_ROBOPTIM")
  PKG_CONFIG_APPEND_CFLAGS (-DBUILDING_ROBOPTIM)
ENDIF()

PKG_CONFIG_APPEND_LIBS(roboptim-core-manifold)

PKG_CONFIG_USE_DEPENDENCY(roboptim-core-manifold roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(roboptim-core-manifold manifolds)

INSTALL(TARGETS roboptim-core-manifold DESTINATION ${CMAKE_INSTALL_LIBDIR})
