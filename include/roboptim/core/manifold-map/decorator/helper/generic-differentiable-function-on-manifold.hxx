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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_GENERIC_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX
#define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_GENERIC_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX

namespace roboptim
{
  template<typename U>
  GenericDifferentiableFunctionOnManifold<U>::GenericDifferentiableFunctionOnManifold(const roboptim::GenericDifferentiableFunction<U>* wrappedFunction, const mnf::Manifold* manifold)
    : wrappedFunction_(wrappedFunction),
      instManifold_(manifold)
  {}

  template<typename U>
  const roboptim::GenericDifferentiableFunction<U>* GenericDifferentiableFunctionOnManifold<U>::getWrappedFunction() const
  {
    return this->wrappedFunction_;
  }

  template<typename U>
  const mnf::Manifold* GenericDifferentiableFunctionOnManifold<U>::getManifold() const
  {
    return this->instManifold_;
  }

}

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_GENERIC_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX
