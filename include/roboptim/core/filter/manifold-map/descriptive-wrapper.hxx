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
        throw std::runtime_error ("Representation dims mismatch on Elementary manifold "
                                  + manifold_->name());
    }
    else
    {
      long size = 0;
      size_t manifolds = manifold_->numberOfSubmanifolds();
      for (size_t i = 0; i < manifolds; ++i)
        size += ((*manifold_) (i)).representationDim();
      if (fct_->inputSize() != size)
        throw std::runtime_error ("Invalid Manifold "
                                  + manifold_->name());
      if (fct_->inputSize() != size)
      {
        std::cerr << "size == " << size << "\n";
        std::cerr << "!= " << fct_->inputSize();
        throw std::runtime_error ("Representation dims mismatch on Composed manifold "
                                  + manifold_->name());
      }
    }
  }

  template <typename U>
  DescriptiveWrapper<U>::~DescriptiveWrapper()
  {
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
