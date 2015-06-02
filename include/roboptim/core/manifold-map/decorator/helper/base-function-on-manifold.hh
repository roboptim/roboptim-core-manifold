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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_BASE_FUNCTION_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_BASE_FUNCTION_ON_MANIFOLD_HH
# include <vector>
# include <utility>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/function.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
# include <roboptim/core/manifold-map/decorator/dispatcher.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>

namespace roboptim
{
  /// \brief Base class of all FunctionOnManifolds
  ///
  /// Implements everything a FunctionOnManifold on any function type needs.
  /// Users should not use this class, as the FunctionOnManifold automatically
  /// create this one if needed.
  ///
  /// \tparam U function type
  template<typename U>
  class BaseFunctionOnManifold : public detail::AutopromoteTrait<U>::T_type
  {
  public:
    ROBOPTIM_FUNCTION_FWD_TYPEDEFS_ (roboptim::GenericFunction<typename U::traits_t>);
    template<typename V, typename W>
    explicit BaseFunctionOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
      :
      detail::AutopromoteTrait<U>::T_type
      (static_cast<size_type>(problemManifold.representationDim()),
       descWrap.fct().outputSize (),
       (boost::format ("%1%")
	% descWrap.fct().getName ()).str ()),
      fct_ (&descWrap.fct()),
      manifold_ (&descWrap.manifold())
    {
      computeMapping(descWrap,
                     problemManifold,
                     functionManifold,
                     restrictedManifolds,
                     restrictions);
    }

  protected:
    void impl_compute (result_ref result, const_argument_ref x)
      const;

    template <typename V, typename W>
    void computeMapping(DescriptiveWrapper<V, W>& descWrap,
			const mnf::Manifold& problemManifold,
			const mnf::Manifold& functionManifold,
			std::vector<const mnf::Manifold*> restrictedManifolds,
			std::vector<std::pair<long, long>> restrictions);

    /// \brief map the input to the restricted problem
    ///
    /// \param argument the argument to map
    void mapArgument(const_argument_ref argument)
      const;

    ~BaseFunctionOnManifold ();

    friend Dispatcher<U, typename U::traits_t>;

    /// \brief the function.
    const U* fct_;
    /// \brief the problem manifold.
    const mnf::Manifold* manifold_;

    /// \brief array representing the restricted mapping
    size_t* mappingFromFunction_;
    /// \brief size of the array
    long mappingFromFunctionSize_;

    /// \brief array representing the restricted mapping for the tangent problem
    size_t* tangentMappingFromFunction_;
    /// \brief size of the array
    long tangentMappingFromFunctionSize_;

    /// \brief new input mapped to the restricted problem
    mutable vector_t mappedInput_;
  };

}

# include <roboptim/core/manifold-map/decorator/helper/base-function-on-manifold.hxx>

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_BASE_FUNCTION_ON_MANIFOLD_HH
