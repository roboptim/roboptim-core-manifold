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
  /// \addtogroup roboptim_decorator
  /// @{

  /// \brief Binds a DescriptiveWrapper to a instance of a submanifold.
  /// \tparam U input function type.
  template <typename U>
  class FunctionOnManifold : public detail::AutopromoteTrait<U>::T_type, private boost::noncopyable
  {
  public:
    typedef typename detail::AutopromoteTrait<U>::T_type parentType_t;
    ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (parentType_t);

    typedef boost::shared_ptr<FunctionOnManifold> FunctionOnManifoldShPtr_t;

    /// \brief Create a mapping from the problem's manifold to the function's.
    /// \param fct input function.
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    /// \param restrictedManifolds a list of elementary Manifolds to be restricted to a part of themselves
    /// \param restrictions the restrictions applying to the selected manifolds, represented as (startingIndex, size). If a single one is given, it will apply to all restricted manifolds.

    // Empty constructor, should be initialized
    // later via a manual call to computeMapping
    explicit FunctionOnManifold();

    template <typename V, typename W>
    explicit FunctionOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const pgs::Manifold& problemManifold,
     const pgs::Manifold& functionManifold,
     std::vector<const pgs::Manifold*> restrictedManifolds,
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

    template<typename V, typename W>
    explicit FunctionOnManifold (DescriptiveWrapper<V, W>& fct,
				 const pgs::Manifold& problemManifold,
				 const pgs::Manifold& functionManifold)
      : FunctionOnManifold(fct, problemManifold, functionManifold,
			   std::vector<const pgs::Manifold*>(), std::vector<std::pair<long, long>>())
    {
    }

    template <typename V, typename W>
    void computeMapping(DescriptiveWrapper<V, W>& descWrap,
			const pgs::Manifold& problemManifold,
			const pgs::Manifold& functionManifold,
			std::vector<const pgs::Manifold*> restrictedManifolds,
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

    void manifold_jacobian (pgs::RefMat jacobian,
			    const_argument_ref arg)
      const;

    std::ostream& print_(std::ostream& o);
  private:
  public:
    const U* fct_;
    pgs::Manifold* manifold_;

    size_t* mappingFromFunction_;
    long mappingFromFunctionSize_;

    size_t* tangentMappingFromFunction_;
    long tangentMappingFromFunctionSize_;

    mutable vector_t mappedInput_;
    mutable derivative_t mappedGradient_;
    mutable jacobian_t mappedJacobian_;

    // Dirty copy ONLY
    mutable Eigen::MatrixXd tangentMappedJacobian;

    void mapArgument(const_argument_ref argument)
      const;

    void unmapGradient(gradient_ref gradient, Eigen::VectorXd& mappedGradient)
      const;

    void unmapTangentJacobian(pgs::RefMat jacobian)
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
