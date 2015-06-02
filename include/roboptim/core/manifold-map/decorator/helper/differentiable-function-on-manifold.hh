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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
# include <vector>
# include <utility>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
# include <roboptim/core/manifold-map/decorator/dispatcher.hh>

# include <manifolds/Manifold.h>

namespace roboptim
{
  /// \brief Base class of FunctionOnManifolds on DifferentiableFunction
  ///
  /// Implements everything a FunctionOnManifold on a DifferentiableFunction
  /// needs, especially the gradient and jacobian support.
  /// Users should not use this class, as the FunctionOnManifold automatically
  /// create this one if needed.
  ///
  /// \tparam U function type
  template <typename U>
  class DifferentiableFunctionOnManifold : public BaseFunctionOnManifold<U>, public GenericDifferentiableFunctionOnManifold<typename U::traits_t>
  {
  public:
    ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (roboptim::GenericDifferentiableFunction<typename U::traits_t>);
  protected:
    DifferentiableFunctionOnManifold() = delete;

  public:
    template<typename V, typename W>
    explicit DifferentiableFunctionOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
      : BaseFunctionOnManifold<U>(descWrap, problemManifold, functionManifold, restrictedManifolds, restrictions),
      GenericDifferentiableFunctionOnManifold<typename U::traits_t> (&descWrap.fct(), &functionManifold)
    {
      this->mappedGradient_ = gradient_t(static_cast<int> (this->mappingFromFunctionSize_));
      this->mappedJacobian_ = jacobian_t(descWrap.fct().outputSize(), static_cast<int> (this->mappingFromFunctionSize_));
      this->tangentMappedJacobian = Eigen::MatrixXd::Zero(descWrap.fct().outputSize(), this->tangentMappingFromFunctionSize_);
      this->mappedGradient_.setZero();
      this->mappedJacobian_.setZero();
    }

    /// \brief apply the jacobian on the manifold's tangent space
    void manifold_jacobian (mnf::RefMat jacobian,
                            const_argument_ref arg)
      const;

    friend Dispatcher<U, typename U::traits_t>;

  protected:
    void impl_gradient (gradient_ref gradient,
                        const_argument_ref argument,
                        size_type functionId = 0)
      const;
    void impl_jacobian (jacobian_ref jacobian,
                        const_argument_ref arg)
      const;

    /// \brief gets the gradient from the restricted problem
    ///
    /// \param gradient the output gradient
    void unmapGradient(gradient_ref) const;

    /// \brief unmap the jacobian from the restricted problem
    ///
    /// \param jacobian the jacobian to unmap
    void unmapTangentJacobian(mnf::RefMat jacobian)
      const;

    /// \brief new gradient mapped to the restricted problem
    mutable gradient_t mappedGradient_;
    /// \brief new jacobian mapped to the restricted problem
    mutable jacobian_t mappedJacobian_;

    // Dirty copy ONLY
    mutable Eigen::MatrixXd tangentMappedJacobian;
  };

}

# include <roboptim/core/manifold-map/decorator/helper/differentiable-function-on-manifold.hxx>

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
