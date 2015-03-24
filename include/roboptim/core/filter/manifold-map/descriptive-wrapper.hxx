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
    if (manifold_->isElementary())
    {
      if (fct_->inputSize() != manifold_->representationDim())
      {
        std::stringstream* error = new std::stringstream;
        (*error) << "Representation dims mismatch on Elementary manifold "
                 << manifold_->name() << " using function "
                 << fct->getName() << ". Expected dimension :"
                 << manifold_->representationDim() << ", actual one is "
                 << fct_->inputSize();
        throw std::runtime_error (error->str());
      }
    }
    else
    {
      long size = 0;
      size_t manifolds = manifold_->numberOfSubmanifolds();
      for (size_t i = 0; i < manifolds; ++i)
        size += ((*manifold_) (i)).representationDim();
      if (manifold_->representationDim() != size)
      {
        std::stringstream* error = new std::stringstream;
        (*error) << "Invalid Composed Manifold " << manifold_->name()
                 << ". Size equals " << manifold_->representationDim()
                 << " whereas the sum of submanifold's equals " << size;
        throw std::runtime_error (error->str());
      }
      if (fct_->inputSize() != size)
      {
        std::stringstream* error = new std::stringstream;
        (*error) << "Representation dims mismatch on Composed manifold "
                 << manifold_->name() << " using function "
                 << fct->getName() << ". Expected dimension :" << size
                 << ", actual one is " << fct_->inputSize();
        throw std::runtime_error (error->str());
      }
    }
  }

  template <typename U>
  DescriptiveWrapper<U>::~DescriptiveWrapper()
  {
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
