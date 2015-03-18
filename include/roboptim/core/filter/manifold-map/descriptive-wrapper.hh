#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
# include <vector>
# include <boost/shared_ptr.hpp>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/differentiable-function.hh>

# include <manifolds/Manifold.h>


namespace roboptim
{
  /// \addtogroup roboptim_filter
  /// @{

  /// \brief apply a given roboptim function to a given dimension space
  /// \tparam U input function type.
  template <typename U>
  class DescriptiveWrapper : public detail::AutopromoteTrait<U>::T_type
  {
  public:
    typedef typename detail::AutopromoteTrait<U>::T_type parentType_t;
    ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (parentType_t);

    typedef boost::shared_ptr<DescriptiveWrapper> DescriptiveWrapperShPtr_t;

    /// \brief set the functions parameters regarding the given manifold.
    /// \param fct input function.
    /// \param functionManifold the manifold describing the function's input vector.
    explicit DescriptiveWrapper (boost::shared_ptr<U> fct,
			  pgs::Manifold& functionManifold);
    ~DescriptiveWrapper ();

    const boost::shared_ptr<U>& fct () const
    {
      return fct_;
    }

    U& fct ()
    {
      return fct_;
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
    boost::shared_ptr<U>  fct_;
    vector_t              x_;
    gradient_t            gradient_;
    jacobian_t            jacobian_;
  };

  template <typename U>
  boost::shared_ptr<DescriptiveWrapper<U> >
  descriptivewrapper (boost::shared_ptr<U> fct,
	  pgs::Manifold& manifold)
  {
    return boost::make_shared<DescriptiveWrapper<U> > (fct, manifold);
  }

} // end of namespace roboptim.

# include <roboptim/core/filter/manifold-map/descriptive-wrapper.hxx>
#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
