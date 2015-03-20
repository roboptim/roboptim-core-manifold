#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
# include <boost/format.hpp>
# include <manifolds/Manifold.h>

namespace roboptim
{

  template <typename U>
  DescriptiveWrapper<U>::DescriptiveWrapper
  (boost::shared_ptr<U> fct,
   pgs::Manifold& functionManifold)
    : detail::AutopromoteTrait<U>::T_type
      (functionManifold.representationDim(),
       fct->outputSize (),
       (boost::format ("%1%")
	% fct->getName ()).str ()),
      fct_ (fct),
      manifold_ (&functionManifold)
  {
    if (manifold_->isElementary())
    {
      if (fct_->inputSize() != manifold_->representationDim())
        throw std::runtime_error ("Representation dims mismatch");
    }
    else
    {
      long size = 0;
      size_t manifolds = manifold_->numberOfSubmanifolds();
      for (size_t i = 0; i < manifolds; ++i)
        size += ((*manifold_) (i)).representationDim();
      if (fct_->inputSize() != size)
      {
        std::cerr << "size == " << size << "\n";
        std::cerr << "!= " << fct_->inputSize();
        throw std::runtime_error ("Representation dims mismatch");
      }
    }
  }

  template <typename U>
  DescriptiveWrapper<U>::~DescriptiveWrapper()
  {
  }

  template <typename U>
  void
  DescriptiveWrapper<U>::impl_compute
  (result_ref result, const_argument_ref x)
    const
  {
    fct_->operator () (result, x);
  }

  template <typename U>
  void
  DescriptiveWrapper<U>::impl_gradient
  (gradient_ref gradient, const_argument_ref arg, size_type functionId)
    const
  {
    fct_->gradient(gradient, arg, functionId);
  }

  template <typename U>
  void
  DescriptiveWrapper<U>::impl_jacobian
  (jacobian_ref jacobian, const_argument_ref arg)
    const
  {
    fct_->jacobian(jacobian, arg);
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
