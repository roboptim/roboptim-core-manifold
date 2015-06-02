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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_TWICE_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_TWICE_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
# include <vector>
# include <utility>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/twice-differentiable-function.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>

# include <manifolds/Manifold.h>

namespace roboptim
{
  /// \brief Base class of FunctionOnManifolds on TwiceDifferentiableFunction
  ///
  /// Implements everything a FunctionOnManifold on a
  /// TwiceDifferentiableFunction needs, that is the hessian support.
  /// Users should not use this class, as the FunctionOnManifold automatically
  /// create this one if needed.
  ///
  /// \tparam U function type
  template <typename U>
  class TwiceDifferentiableFunctionOnManifold : public DifferentiableFunctionOnManifold<U>
  {
  public:
    ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (roboptim::GenericTwiceDifferentiableFunction<typename U::traits_t>);
  protected:
    TwiceDifferentiableFunctionOnManifold() = delete;

  public:
    template <typename V, typename W>
    explicit TwiceDifferentiableFunctionOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
      : DifferentiableFunctionOnManifold<U>(descWrap, problemManifold, functionManifold, restrictedManifolds, restrictions)
    {
      this->mappedHessian_ = hessian_t(descWrap.fct().inputSize(), static_cast<int> (this->mappingFromFunctionSize_));
      this->mappedHessian_.setZero();
    }

  protected:
    void impl_hessian(hessian_ref, const_argument_ref, size_type) const;

    /// \brief gets the hessian from the restricted problem
    ///
    /// \param hessian the output hessian
    void unmapHessian(hessian_ref hessian) const;

    friend Dispatcher<U, typename U::traits_t>;

    /// \brief new hessian mapped to the restricted problem
    mutable hessian_t mappedHessian_;
  };
}

# include <roboptim/core/manifold-map/decorator/helper/twice-differentiable-function-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_TWICE_DIFFERENTIABLE_FUNCTION_ON_MANIFOLD_HH
