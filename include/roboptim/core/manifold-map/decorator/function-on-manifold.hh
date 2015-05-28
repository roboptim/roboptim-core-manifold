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
# include <type_traits>

# include <boost/shared_ptr.hpp>
# include <boost/make_shared.hpp>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/twice-differentiable-function.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
# include <roboptim/core/manifold-map/decorator/dispatcher.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>
# include <manifolds/S2.h>

# include <boost/noncopyable.hpp>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  template<typename U>
  class GenericDifferentiableFunctionOnManifold
  {
  protected:
    GenericDifferentiableFunctionOnManifold() = delete;

    const roboptim::GenericDifferentiableFunction<U>* wrappedFunction_;
    const mnf::Manifold* manifold_;
  public:
    GenericDifferentiableFunctionOnManifold(const roboptim::GenericDifferentiableFunction<U>* wrappedFunction, const mnf::Manifold* manifold);

    const roboptim::GenericDifferentiableFunction<U>* getWrappedFunction() const;
    const mnf::Manifold* getManifold() const;
  };

  // TODO: Move to a separate file
  template<typename U>
  class BecauseUsersCanDoShit
  {
  public:
    ROBOPTIM_FUNCTION_FWD_TYPEDEFS_ (roboptim::GenericFunction<typename U::traits_t>);
    template<typename V, typename W>
    explicit BecauseUsersCanDoShit
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
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

  template <typename U>
  class DifferentiableFunctionOnManifold : public BecauseUsersCanDoShit<U>, public GenericDifferentiableFunctionOnManifold<typename U::traits_t>
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
      : GenericDifferentiableFunctionOnManifold<typename U::traits_t> (&descWrap.fct(), &functionManifold)
    {
      computeMapping(descWrap,
                     problemManifold,
                     functionManifold,
                     restrictedManifolds,
                     restrictions);
    }

  protected:
    template <typename V, typename W>
    void computeMapping(DescriptiveWrapper<V, W>& descWrap,
			const mnf::Manifold& problemManifold,
			const mnf::Manifold& functionManifold,
			std::vector<const mnf::Manifold*> restrictedManifolds,
			std::vector<std::pair<long, long>> restrictions);

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

    /// \brief gets the gradient from the restricted problem
    ///
    /// \param instance the FunctionOnManifold
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

  template <typename U>
  class TwiceDifferentiableFunctionOnManifold : public DifferentiableFunctionOnManifold<U>
  {
  public:
    ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (roboptim::GenericTwiceDifferentiableFunction<typename U::traits_t>);
  protected:
    TwiceDifferentiableFunctionOnManifold() = delete;

  public:
    template <typename V, typename W>
    void computeMapping(DescriptiveWrapper<V, W>& descWrap,
			const mnf::Manifold& problemManifold,
			const mnf::Manifold& functionManifold,
			std::vector<const mnf::Manifold*> restrictedManifolds,
			std::vector<std::pair<long, long>> restrictions);

    template <typename V, typename W>
    explicit TwiceDifferentiableFunctionOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
    {
      computeMapping(descWrap,
                     problemManifold,
                     functionManifold,
                     restrictedManifolds,
                     restrictions);
    }

  protected:
    void impl_hessian(hessian_ref, const_argument_ref, size_type) const;
    void unmapHessian(hessian_ref hessian) const;

    mutable hessian_t mappedHessian_;
  };

  /// \brief Maps a DescriptiveWrapper to a instance of a submanifold.
  ///
  /// \tparam U input roboptim function type.
  // TODO: deactivate gradient and jacobian when not inheriting from DifferentiableFunction
  template <typename U>
  class FunctionOnManifold :
    public detail::AutopromoteTrait<U>::T_type,
    private boost::noncopyable,
    public std::conditional<std::is_base_of<roboptim::GenericTwiceDifferentiableFunction<typename U::traits_t>,
					    U>::value,
			    TwiceDifferentiableFunctionOnManifold<U>,
			    typename std::conditional<std::is_base_of<roboptim::GenericDifferentiableFunction<typename U::traits_t>,
							     U>::value,
					     DifferentiableFunctionOnManifold<U>,
					     BecauseUsersCanDoShit<U> >::type >::type
  {
  public:
    ROBOPTIM_FUNCTION_FWD_TYPEDEFS_ (U);
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
	(static_cast<size_type>(problemManifold.representationDim()),
	 descWrap.fct().outputSize (),
	 (boost::format ("%1%")
	  % descWrap.fct().getName ()).str ()),
	std::conditional<std::is_base_of<roboptim::GenericTwiceDifferentiableFunction<typename U::traits_t>,
					 U>::value,
			 TwiceDifferentiableFunctionOnManifold<U>,
			 typename std::conditional<std::is_base_of<roboptim::GenericDifferentiableFunction<typename U::traits_t>,
							  U>::value,
					  DifferentiableFunctionOnManifold<U>,
					  BecauseUsersCanDoShit<U> >::type >::type(descWrap, problemManifold, functionManifold, restrictedManifolds, restrictions)
    {
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

    /// \brief Traits type.
    typedef typename parent_t::traits_t traits_t;

    ~FunctionOnManifold ();

    std::ostream& print_(std::ostream& o);
  private:
  public:
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
