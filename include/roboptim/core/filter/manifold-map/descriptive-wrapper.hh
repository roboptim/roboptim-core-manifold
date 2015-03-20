#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
# include <vector>
# include <ostream>
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

    /// \brief map the function on the given manifold.
    /// \param fct input function.
    /// \param functionManifold the manifold describing the function's input vector.
    explicit DescriptiveWrapper (boost::shared_ptr<U> fct,
			  pgs::Manifold& functionManifold);
    ~DescriptiveWrapper ();

    const pgs::Manifold& manifold () const
    {
      return *manifold_;
    }

    pgs::Manifold& manifold ()
    {
      return *manifold_;
    }

    const boost::shared_ptr<U>& fct () const
    {
      return fct_;
    }

    U& fct ()
    {
      return *fct_;
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
    pgs::Manifold* manifold_;
  };

  template <typename U>
  boost::shared_ptr<DescriptiveWrapper<U> >
  descriptivewrapper (boost::shared_ptr<U> fct,
	  pgs::Manifold& manifold)
  {
    return boost::make_shared<DescriptiveWrapper<U> > (fct, manifold);
  }

  template <typename U>
  std::ostream&
  operator<<(std::ostream& o, DescriptiveWrapper<U>& descWrap)
  {
    o <<
    "Displaying info about function in DescriptiveWrapper :" << "\n" <<
    "the function is " << descWrap.fct().getName() << "\n" <<
    "the function input size is " << descWrap.fct().inputSize() << "\n" <<
    "the function output size is " << descWrap.fct().outputSize() << "\n" <<
    "Displaying info about the function manifold in DescriptiveWrapper :" << "\n" <<
    "the manifold is " << descWrap.manifold().name() << "\n" <<
    "Is the manifold elementary ? : " <<
    (descWrap.manifold().isElementary() ? "yes" : "no" ) << "\n" <<
    "manifold display :" << "\n";
    descWrap.manifold().display();
    o << "\n" <<
    "the manifold dimension is " << descWrap.manifold().representationDim() << "\n";

    return o;
  }
} // end of namespace roboptim.

# include <roboptim/core/filter/manifold-map/descriptive-wrapper.hxx>
#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
