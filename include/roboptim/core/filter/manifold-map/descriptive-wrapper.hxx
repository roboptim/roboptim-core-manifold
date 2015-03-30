#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
# include <boost/format.hpp>
# include <manifolds/Manifold.h>

namespace roboptim
{

  template <typename U, typename V>
  template<class ... Types>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (Types ... args)
  {
    this->fct_ = new U(args...);
    V manifoldGenerator;
    this->manifold_ = manifoldGenerator.getManifold(this->function_);
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (boost::shared_ptr<U>& fct, const V& manifold)
  : fct_ (fct),
    manifold_ (manifold)
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

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
