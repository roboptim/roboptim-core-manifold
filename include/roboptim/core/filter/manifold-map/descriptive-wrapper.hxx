#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
# include <boost/format.hpp>
# include <manifolds/Manifold.h>

namespace roboptim
{

  template <typename U>
  DescriptiveWrapper<U>::DescriptiveWrapper
  (boost::shared_ptr<U> fct,
   const pgs::Manifold& functionManifold)
    : fct_ (fct),
      manifold_ (&functionManifold)
  {
    long size = manifold_->representationDim();
    if (fct_->inputSize() != size)
    {
      std::stringstream* error = new std::stringstream;
      (*error) << "Representation dims mismatch on manifold "
               << manifold_->name() << " using function "
               << fct->getName() << ". Expected dimension :" << size
               << ", actual one is " << fct_->inputSize();
      throw std::runtime_error (error->str());
    }
  }

  template <typename U>
  DescriptiveWrapper<U>::~DescriptiveWrapper()
  {
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
