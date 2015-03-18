#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HH
# include <vector>
# include <iostream>
# include <boost/shared_ptr.hpp>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/differentiable-function.hh>

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
    explicit InstanceWrapper (boost::shared_ptr<U> fct,
			  pgs::Manifold& problemManifold,
			  pgs::Manifold& functionManifold);
    ~InstanceWrapper ();

    const boost::shared_ptr<U>& origin () const
    {
      return origin_;
    }

    boost::shared_ptr<U>& origin ()
    {
      return origin_;
    }

    void impl_compute (result_ref result, const_argument_ref x)
      const;

    void impl_gradient (gradient_ref gradient,
			const_argument_ref argument,
			size_type functionId = 0)
      const;
    void impl_jacobian (jacobian_ref jacobian,
			const_argument_ref arg)
      const;

    std::ostream& print(std::ostream& o);
  private:
    boost::shared_ptr<U> origin_;

    int* mappingFromProblem_;
    size_t mappingFromProblemSize_;
    int* mappingFromFunction_;
    size_t mappingFromFunctionSize_;
  };

  template <typename U>
  boost::shared_ptr<InstanceWrapper<U> >
  scalar (boost::shared_ptr<U> origin,
	  typename InstanceWrapper<U>::size_type start = 0,
	  typename InstanceWrapper<U>::size_type size = 1)
  {
    return boost::make_shared<InstanceWrapper<U> > (origin, start, size);
  }

  template <typename U>
  boost::shared_ptr<InstanceWrapper<U> >
  operator* (typename InstanceWrapper<U>::value_type scalar,
	     boost::shared_ptr<U> origin)
  {
    return boost::make_shared<InstanceWrapper<U> > (origin, scalar);
  }

  template <typename U>
  boost::shared_ptr<InstanceWrapper<U> >
  operator* (boost::shared_ptr<U> origin,
	     typename InstanceWrapper<U>::value_type scalar)
  {
    return boost::make_shared<InstanceWrapper<U> > (origin, scalar);
  }

  template <typename U>
  boost::shared_ptr<U>
  operator+ (boost::shared_ptr<U> origin)
  {
    return origin;
  }

  template <typename U>
  boost::shared_ptr<InstanceWrapper<U> >
  operator- (boost::shared_ptr<U> origin)
  {
    return boost::make_shared<InstanceWrapper<U> > (origin, -1.);
  }

  template <typename U>
  std::ostream&
  operator<<(std::ostream& o, InstanceWrapper<U>& instWrap)
  {
    return instWrap.print(o);
  }

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/filter/manifold-map/instance-wrapper.hxx>
#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HH
