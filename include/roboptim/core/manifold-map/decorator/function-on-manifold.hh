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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
# include <vector>
# include <iostream>
# include <utility>
# include <boost/shared_ptr.hpp>
# include <boost/make_shared.hpp>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>
# include <manifolds/S2.h>

# include <boost/noncopyable.hpp>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  /// \brief Maps a DescriptiveWrapper to a instance of a submanifold.
  ///
  /// \tparam U input roboptim function type.
  template <typename U>
  class FunctionOnManifold : public detail::AutopromoteTrait<U>::T_type, private boost::noncopyable
  {
  public:
    typedef typename detail::AutopromoteTrait<U>::T_type parentType_t;
    ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (parentType_t);

    typedef boost::shared_ptr<FunctionOnManifold> FunctionOnManifoldShPtr_t;

    /// \brief Creates the mapping between the manifolds.
    ///
    /// Programmers should ideally not directly use this constructor, amd
    /// rely on the macros defined in the file manifold-map.hh to get the well
    /// defined FunctionOnManifold type, ready to be instantiated.
    ///
    /// \tparam V input function type templating the DescriptiveWrapper.
    /// \tparam W descriptive manifold type templating the DescriptiveWrapper.
    ///
    /// \param descwrap instance of the DescriptiveWrapper
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    /// \param restrictedManifolds a list of elementary Manifolds to be restricted to a part of themselves
    /// \param restrictions the restrictions applying to the selected manifolds, represented as (startingIndex, size). If a single one is given, it will apply to all restricted manifolds.
    template <typename V, typename W>
    explicit FunctionOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
      : detail::AutopromoteTrait<U>::T_type
	(problemManifold.representationDim(),
	 descWrap.fct().outputSize (),
	 (boost::format ("%1%")
	  % descWrap.fct().getName ()).str ()),
	fct_(&descWrap.fct()),
	manifold_(&descWrap.manifold())
    {
      computeMapping(descWrap,
		     problemManifold,
		     functionManifold,
		     restrictedManifolds,
		     restrictions);
    }


    /// \brief Creates the mapping between the manifolds without restrictions.
    ///
    /// Programmers should ideally not directly use this constructor, amd
    /// rely on the macros defined in the file manifold-map.hh to get the well
    /// defined FunctionOnManifold type, ready to be instantiated.
    ///
    /// \tparam V input function type templating the DescriptiveWrapper.
    /// \tparam W descriptive manifold type templating the DescriptiveWrapper.
    ///
    /// \param fct instance of the DescriptiveWrapper
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    template<typename V, typename W>
    explicit FunctionOnManifold (DescriptiveWrapper<V, W>& fct,
				 const mnf::Manifold& problemManifold,
				 const mnf::Manifold& functionManifold)
      : FunctionOnManifold(fct, problemManifold, functionManifold,
			   std::vector<const mnf::Manifold*>(), std::vector<std::pair<long, long>>())
    {
    }

    template <typename V, typename W>
    void computeMapping(DescriptiveWrapper<V, W>& descWrap,
			const mnf::Manifold& problemManifold,
			const mnf::Manifold& functionManifold,
			std::vector<const mnf::Manifold*> restrictedManifolds,
			std::vector<std::pair<long, long>> restrictions);

    /// \brief Traits type.
    typedef typename parent_t::traits_t traits_t;

    /// \brief Derivative type.
    ///
    /// Derivatives are column vectors.
    ROBOPTIM_GENERATE_TRAITS_REFS_T(derivative,traits_t);

    ~FunctionOnManifold ();

    void impl_compute (result_ref result, const_argument_ref x)
      const;

    void impl_gradient (gradient_ref gradient,
			const_argument_ref argument,
			size_type functionId = 0)
      const;
    void impl_jacobian (jacobian_ref jacobian,
			const_argument_ref arg)
      const;

    void manifold_jacobian (mnf::RefMat jacobian,
			    const_argument_ref arg)
      const;

    std::ostream& print_(std::ostream& o);
  private:
  public:
    /// \brief the function.
    const U* fct_;
    /// \brief the problem manifold.
    mnf::Manifold* manifold_;

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
    /// \brief new gradient mapped to the restricted problem
    mutable derivative_t mappedGradient_;
    /// \brief new jacobian mapped to the restricted problem
    mutable jacobian_t mappedJacobian_;

    // Dirty copy ONLY
    mutable Eigen::MatrixXd tangentMappedJacobian;

    /// \brief map the input to the restricted problem
    ///
    /// \param argument the argument to map
    void mapArgument(const_argument_ref argument)
      const;

    /// \brief gets the gradient from the restricted problem
    ///
    /// \param gradient the output gradient
    /// \param mappedGradient the gradient computed to unmap
    void unmapGradient(gradient_ref gradient, const_gradient_ref mappedGradient)
      const;

    /// \brief unmap the jacobian from the restricted problem
    ///
    /// \param jacobian the jacobian to unmap
    void unmapTangentJacobian(mnf::RefMat jacobian)
      const;
  };

  template <typename U>
  std::ostream&
  operator<<(std::ostream& o, FunctionOnManifold<U>& instWrap)
  {
    return instWrap.print_(o);
  }

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/manifold-map/decorator/function-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
