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
  }

  template <typename U>
  void
  DescriptiveWrapper<U>::impl_gradient
  (gradient_ref result, const_argument_ref x, size_type functionId)
    const
  {
  }

  template <typename U>
  void
  DescriptiveWrapper<U>::impl_jacobian
  (jacobian_ref result, const_argument_ref x)
    const
  {
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_DESCRIPTIVE_WRAPPER_HXX
