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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX
#define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX

namespace roboptim
{
  template <typename U>
  void
  DifferentiableFunctionOnManifold<U>::unmapTangentJacobian(mnf::RefMat jacobian)
    const
  {
    for (long i = 0; i < this->tangentMappingFromFunctionSize_; ++i)
      {
	jacobian.col(static_cast<long>(this->tangentMappingFromFunction_[i])) = this->tangentMappedJacobian.col(i);
      }
  }

  template <typename U>
  void
  DifferentiableFunctionOnManifold<U>::unmapGradient(gradient_ref gradient) const
  {
    assert(gradient.cols() == this->inputSize());

    for (int i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	gradient.coeffRef(static_cast<size_type>(this->mappingFromFunction_[i])) = this->mappedGradient_.coeffRef(i);
      }
  }

  template <typename U>
  void
  DifferentiableFunctionOnManifold<U>::impl_gradient (gradient_ref gradient,
						      const_argument_ref argument,
						      size_type functionId)
    const
  {
    this->mapArgument(argument);
    this->fct_->gradient(this->mappedGradient_, this->mappedInput_, functionId);

    gradient.setZero();
    unmapGradient(gradient);
  }

  template <typename U>
  void
  DifferentiableFunctionOnManifold<U>::impl_jacobian (jacobian_ref jacobian,
						      const_argument_ref argument)
    const
  {
    this->mapArgument(argument);
    jacobian.setZero();

    this->fct_->jacobian(this->mappedJacobian_, this->mappedInput_);
    roboptim::Dispatcher<U, typename U::traits_t>::unmapJacobian(this, jacobian);
  }

  template <typename U>
  void
  DifferentiableFunctionOnManifold<U>::manifold_jacobian (mnf::RefMat jacobian,
							  const_argument_ref argument)
    const
  {
    this->mapArgument(argument);
    jacobian.setZero();
    this->mappedJacobian_.setZero();
    this->tangentMappedJacobian.setZero();

    this->fct_->jacobian(this->mappedJacobian_, this->mappedInput_);

    roboptim::Dispatcher<U, typename U::traits_t>::applyDiff(this);

    this->unmapTangentJacobian(jacobian);
  }
}

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HXX
