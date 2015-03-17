#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_HH
# include <vector>
# include <boost/shared_ptr.hpp>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/differentiable-function.hh>

# include <manifolds/Manifold.h>


namespace roboptim
{
  /// \addtogroup roboptim_filter
  /// @{

  /// \brief Expose only a subset of all the variables to the function.
  /// \tparam U input function type.
  template <typename U>
  class ManifoldMap : public detail::AutopromoteTrait<U>::T_type
  {
  public:
    typedef typename detail::AutopromoteTrait<U>::T_type parentType_t;
    ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (parentType_t);

    typedef boost::shared_ptr<ManifoldMap> ManifoldMapShPtr_t;

    /// \brief Create a mapping from the problem's manifold to the function's.
    /// \param fct input function.
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    explicit ManifoldMap (boost::shared_ptr<U> fct,
			  pgs::Manifold& problemManifold,
			  pgs::Manifold& functionManifold);
    ~ManifoldMap ();

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
  private:
    boost::shared_ptr<U> origin_;

    int* mapping_;
    size_t mappingSize_;
  };

  template <typename U>
  boost::shared_ptr<ManifoldMap<U> >
  scalar (boost::shared_ptr<U> origin,
	  typename ManifoldMap<U>::size_type start = 0,
	  typename ManifoldMap<U>::size_type size = 1)
  {
    return boost::make_shared<ManifoldMap<U> > (origin, start, size);
  }

  template <typename U>
  boost::shared_ptr<ManifoldMap<U> >
  operator* (typename ManifoldMap<U>::value_type scalar,
	     boost::shared_ptr<U> origin)
  {
    return boost::make_shared<ManifoldMap<U> > (origin, scalar);
  }

  template <typename U>
  boost::shared_ptr<ManifoldMap<U> >
  operator* (boost::shared_ptr<U> origin,
	     typename ManifoldMap<U>::value_type scalar)
  {
    return boost::make_shared<ManifoldMap<U> > (origin, scalar);
  }

  template <typename U>
  boost::shared_ptr<U>
  operator+ (boost::shared_ptr<U> origin)
  {
    return origin;
  }

  template <typename U>
  boost::shared_ptr<ManifoldMap<U> >
  operator- (boost::shared_ptr<U> origin)
  {
    return boost::make_shared<ManifoldMap<U> > (origin, -1.);
  }

  /// @}

} // end of namespace roboptim.

# include <roboptim/core/filter/manifold-map.hxx>
#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_HH
