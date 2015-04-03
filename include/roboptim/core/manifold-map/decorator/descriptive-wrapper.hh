#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
# include <vector>
# include <ostream>
# include <boost/shared_ptr.hpp>

# include <roboptim/core/differentiable-function.hh>

# include <manifolds/Manifold.h>

namespace roboptim
{
  /// \addtogroup roboptim_decorator
  /// @{

  // FIXME: private constructor MAGIC

  /// \brief apply a given roboptim function to a given dimension space
  /// \tparam U input function type.
  template <typename U, typename V>
  class DescriptiveWrapper
  {
  public:
    typedef boost::shared_ptr<DescriptiveWrapper> DescriptiveWrapperShPtr_t;

    /// \brief map the function on the given manifold.
    /// \param fct input function.
    /// \param functionManifold the manifold describing the function's input vector.
    template<class ... Types>
    explicit DescriptiveWrapper (Types ... args);

    DescriptiveWrapper (boost::shared_ptr<U>& f, const pgs::Manifold& m);

    ~DescriptiveWrapper ();

    const pgs::Manifold& manifold () const
    {
      return *manifold_;
    }

    const boost::shared_ptr<U>& fctPointer () const
    {
      return fct_;
    }

    U& fct ()
    {
      return *fct_;
    }

    const boost::shared_ptr<pgs::Manifold>& manifoldPointer () const
    {
      return manifold_;
    }

  private:

    boost::shared_ptr<U>  fct_;
    boost::shared_ptr<pgs::Manifold> manifold_;
  };

  template <typename U, typename V>
  boost::shared_ptr<DescriptiveWrapper<U, V> >
  descriptivewrapper (boost::shared_ptr<U> fct,
		      const pgs::Manifold& manifold)
  {
    return boost::make_shared<DescriptiveWrapper<U, V> > (fct, manifold);
  }

  template <typename U, typename V>
  std::ostream&
  operator<<(std::ostream& o, DescriptiveWrapper<U, V>& descWrap)
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
    "the manifold dimension is " << descWrap.manifold().representationDim() << "\n";
    return o;
  }

} // end of namespace roboptim.

# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
