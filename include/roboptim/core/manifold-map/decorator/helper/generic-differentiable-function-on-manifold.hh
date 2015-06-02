// Copyright (C) 2015 by Félix Darricau, AIST, CNRS, EPITA
//                       Grégoire Duchemin, AIST, CNRS, EPITA
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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_GENERIC_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_GENERIC_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
# include <roboptim/core/differentiable-function.hh>
# include <manifolds/Manifold.h>

namespace roboptim
{
  /// \brief Interface to manage dynamic casting on FunctionOnManifolds
  ///
  /// \tparam U Eigen matrix type used
  template<typename U>
  class GenericDifferentiableFunctionOnManifold
  {
  protected:
    GenericDifferentiableFunctionOnManifold() = delete;

    const roboptim::GenericDifferentiableFunction<U>* wrappedFunction_;
    const mnf::Manifold* instManifold_;
  public:
    GenericDifferentiableFunctionOnManifold(const roboptim::GenericDifferentiableFunction<U>* wrappedFunction, const mnf::Manifold* manifold);

    const roboptim::GenericDifferentiableFunction<U>* getWrappedFunction() const;
    const mnf::Manifold* getManifold() const;
  };
}

# include <roboptim/core/manifold-map/decorator/helper/generic-differentiable-function-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_GENERIC_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
