#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
# include <vector>
# include <ostream>
# include <boost/shared_ptr.hpp>

# include <roboptim/core/differentiable-function.hh>

# include <manifolds/Manifold.h>

#define BIND_FUNCTION_ON_MANIFOLD(function, manifold) typedef DW<function, manifold> function##_On_##manifold ;
#define DEFAULT_FUNCTION_BINDING(function, manifold) typedef DW<function, manifold> Wrapped##function ;

namespace roboptim
{
  /// \addtogroup roboptim_filter
  /// @{

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

  private:
    DescriptiveWrapper (boost::shared_ptr<U>& f, const Manifold& m);

  public:
    DescriptiveWrapper (boost::shared_ptr<U>& f, const Manifold& m)
    {
      DescriptiveWrapper(f, m);
    }

    ~DescriptiveWrapper ();

    const pgs::Manifold& manifold () const
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

  private:

    boost::shared_ptr<U>  fct_;
    const pgs::Manifold* manifold_;
  };

  template <typename U>
  boost::shared_ptr<DescriptiveWrapper<U> >
  descriptivewrapper (boost::shared_ptr<U> fct,
	  const pgs::Manifold& manifold)
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
    "the manifold dimension is " << descWrap.manifold().representationDim() << "\n";
    return o;
  }

} // end of namespace roboptim.

# include <roboptim/core/filter/manifold-map/descriptive-wrapper.hxx>
#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HH
