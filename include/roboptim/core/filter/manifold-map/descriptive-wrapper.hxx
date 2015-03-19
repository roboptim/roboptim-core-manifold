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
      x_ (fct->inputSize ()),
      gradient_ (fct->inputSize ()),
      jacobian_ (fct->outputSize (),
		 fct->inputSize ())
  {
    if (fct_->inputSize() != functionManifold.representationDim())
      throw std::runtime_error ("Representation dims mismatch");
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
