// Copyright (C) 2015 by Grégoire Duchemin, AIST, CNRS, EPITA
//                       Félix Darricau, AIST, CNRS, EPITA
//
// This file is part of the roboptim.
//
// roboptim is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_DESC_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_DESC_HH



namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  /// \brief Actual descriptive manifold class for non elementary manifold
  ///
  /// Programmer should rely on the DESC_MANIFOLD macro defined in the file
  /// manifold-map.hh to properly define it.
  ///
  /// \tparam Types list of types representing the elementary manifolds
  /// composing the descriptive manifold
  template<template <typename> class ... Types>
  struct ManiDesc
  {
    /// \brief creates the actual manifold from the underlying library, given
    /// a function
    ///
    /// \tparam U the function type
    ///
    /// \param function the function instance.
    template<class U>
    static pgs::Manifold* getManifold(U* function = nullptr);

  };

  /// @}
}

# include "manifold-desc.hxx"

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_DESC_HH
