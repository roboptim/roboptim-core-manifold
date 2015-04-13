#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HXX
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HXX
# include <boost/format.hpp>
# include <manifolds/Manifold.h>
# include "manifold-desc.hh"

namespace roboptim
{

  template <typename U, typename V>
  template<class ... Types>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (Types ... args)
  {
    this->fct_ = new U(args...);
    this->manifold_ = V::getManifold(this->fct_);
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (U* fct, pgs::Manifold& manifold)
  : fct_ (fct),
    manifold_ (&manifold)
  {
    long size = manifold_->representationDim();
    if (fct_->inputSize() != size)
    {
      std::stringstream* error = new std::stringstream;
      (*error) << "Representation dims mismatch on manifold "
               << manifold_->name() << " using function "
               << fct_->getName() << ". Expected dimension :" << size
               << ", actual one is " << fct_->inputSize();
      throw std::runtime_error (error->str());
    }
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::~DescriptiveWrapper()
  {
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HXX
