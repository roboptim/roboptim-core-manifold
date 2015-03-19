#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HH
# include <vector>
# include <iostream>
# include <boost/shared_ptr.hpp>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/filter/manifold-map/descriptive-wrapper.hh>

# include <manifolds/Manifold.h>


namespace roboptim
{
  /// \addtogroup roboptim_filter
  /// @{

  /// \brief Binds a DescriptiveWrapper to a instance of a submanifold.
  /// \tparam U input function type.
  template <typename U>
  class InstanceWrapper : public detail::AutopromoteTrait<U>::T_type
  {
  public:
    typedef typename detail::AutopromoteTrait<U>::T_type parentType_t;
    ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (parentType_t);

    typedef boost::shared_ptr<InstanceWrapper> InstanceWrapperShPtr_t;

    /// \brief Create a mapping from the problem's manifold to the function's.
    /// \param fct input function.
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    explicit InstanceWrapper (boost::shared_ptr<DescriptiveWrapper<U>> fct,
			      pgs::Manifold& problemManifold,
			      pgs::Manifold& functionManifold);
    ~InstanceWrapper ();

    void impl_compute (result_ref result, const_argument_ref x)
      const;

    void impl_gradient (gradient_ref gradient,
			const_argument_ref argument,
			size_type functionId = 0)
      const;
    void impl_jacobian (jacobian_ref jacobian,
			const_argument_ref arg)
      const;

    std::ostream& print_(std::ostream& o);
  private:
  public:
    boost::shared_ptr<DescriptiveWrapper<U>> descWrap_;

    size_t* mappingFromProblem_;
    long mappingFromProblemSize_;
    size_t* mappingFromFunction_;
    long mappingFromFunctionSize_;

    mutable vector_t mappedInput_;
    mutable gradient_t mappedGradient_;
    mutable jacobian_t mappedJacobian_;

    void mapArgument(const_argument_ref argument)
      const;

    void unmapGradient(gradient_ref gradient)
      const;
  };

  template <typename U>
  boost::shared_ptr<InstanceWrapper<U> >
  scalar (boost::shared_ptr<U> origin,
	  typename InstanceWrapper<U>::size_type start = 0,
	  typename InstanceWrapper<U>::size_type size = 1)
  {
    return boost::make_shared<InstanceWrapper<U> > (origin->fct(), start, size);
  }

  template <typename U>
  boost::shared_ptr<InstanceWrapper<U> >
  operator* (typename InstanceWrapper<U>::value_type scalar,
	     boost::shared_ptr<U> origin)
  {
    return boost::make_shared<InstanceWrapper<U> > (origin->fct(), scalar);
  }

  template <typename U>
  boost::shared_ptr<InstanceWrapper<U> >
  operator* (boost::shared_ptr<U> origin,
	     typename InstanceWrapper<U>::value_type scalar)
  {
    return boost::make_shared<InstanceWrapper<U> > (origin->fct(), scalar);
  }

  template <typename U>
  boost::shared_ptr<U>
  operator+ (boost::shared_ptr<U> origin)
  {
    return origin->fct();
  }

  template <typename U>
  boost::shared_ptr<InstanceWrapper<U> >
  operator- (boost::shared_ptr<U> origin)
  {
    return boost::make_shared<InstanceWrapper<U> > (origin->fct(), -1.);
  }

  template <typename U>
  std::ostream&
  operator<<(std::ostream& o, InstanceWrapper<U>& instWrap)
  {
    return instWrap.print_(o);
  }

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/filter/manifold-map/instance-wrapper.hxx>
#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HH
