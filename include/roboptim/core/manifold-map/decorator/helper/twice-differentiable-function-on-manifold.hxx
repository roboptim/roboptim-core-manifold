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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_TWICE_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX
#define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_TWICE_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX

namespace roboptim
{
  template<typename U>
  void
  TwiceDifferentiableFunctionOnManifold<U>::unmapHessian(hessian_ref hessian) const
  {
    assert (hessian.rows() == this->inputSize());
    assert (hessian.cols() == this->inputSize());

    for (int i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	for (int j = 0; j < this->inputSize(); ++j)
	  hessian.coeffRef(i,j) = this->mappedHessian_.coeffRef(i,j);
      }
  }

  template<typename U>
  void
  TwiceDifferentiableFunctionOnManifold<U>::impl_hessian(hessian_ref hessian, const_argument_ref argument, size_type functionId) const
  {
    this->mapArgument(argument);
    hessian.setZero();

    this->fct_->hessian(this->mappedHessian_, this->mappedInput_, functionId);
    unmapHessian(hessian);
  }
}

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_TWICE_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX
