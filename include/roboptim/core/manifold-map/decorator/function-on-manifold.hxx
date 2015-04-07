#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HXX
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HXX
# include <vector>
# include <queue>
# include <utility>
# include <iostream>
# include <typeinfo>

# include <boost/format.hpp>
# include <boost/mpl/assert.hpp>

# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>
# include <manifolds/S2.h>
namespace roboptim
{

  template <typename U>
  FunctionOnManifold<U>::~FunctionOnManifold()

  {
    delete [] this->mappingFromFunction_;
  }

  template <typename U>
  void
  FunctionOnManifold<U>::mapArgument (const_argument_ref argument)
    const
  {
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	this->mappedInput_(i) = argument(static_cast<long>(this->mappingFromFunction_[i]));
      }
  }

  template <typename U>
  void
  FunctionOnManifold<U>::unmapGradient(gradient_ref gradient, Eigen::VectorXd& mappedGradient)
    const
  {
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	gradient(static_cast<long>(this->mappingFromFunction_[i])) = mappedGradient(i);
      }
  }

  template <typename U>
  void
  FunctionOnManifold<U>::unmapTangentJacobian(pgs::RefMat jacobian)
  const
  {
    for (long i = 0; i < this->tangentMappingFromFunctionSize_; ++i)
    {
      jacobian.col(static_cast<long>(this->tangentMappingFromFunction_[i])) = this->tangentMappedJacobian.col(i);
    }
  }

  template <typename U>
  void
  FunctionOnManifold<U>::impl_compute
  (result_ref result, const_argument_ref x)
    const
  {
    this->mapArgument(x);
    (*this->fct_)(result, this->mappedInput_);
  }

  template <typename U>
  void
  FunctionOnManifold<U>::impl_gradient (gradient_ref gradient,
			 const_argument_ref argument,
			 size_type functionId)
    const
  {
    this->mapArgument(argument);
    this->fct_->gradient(this->mappedGradient_, this->mappedInput_, functionId);

    gradient.setZero();
    this->unmapGradient(gradient, this->mappedGradient_);
  }

  template <typename U>
  void
  FunctionOnManifold<U>::impl_jacobian (jacobian_ref jacobian,
			 const_argument_ref argument)
    const
  {
    this->mapArgument(argument);
    jacobian.setZero();

    for (long j = 0; j < jacobian.rows(); ++j)
      {
	this->fct_->gradient(this->mappedGradient_, this->mappedInput_, j);
	this->unmapGradient(jacobian.row(j), this->mappedGradient_);
      }
  }

  template <typename U>
  void
  FunctionOnManifold<U>::manifold_jacobian (pgs::RefMat jacobian,
                                            const_argument_ref argument)\
  const
  {
    this->mapArgument(argument);
    jacobian.setZero();
    this->mappedJacobian_.setZero();
    this->tangentMappedJacobian.setZero();

    this->fct_->jacobian(this->mappedJacobian_, this->mappedInput_);

    this->manifold_->applyDiffRetractation(this->tangentMappedJacobian, this->mappedJacobian_, this->mappedInput_);

    this->unmapTangentJacobian(jacobian);
  }

  template <typename U>
  std::ostream&
  FunctionOnManifold<U>::print_ (std::ostream& o)
  {
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	o << (i>0?", ":"") << this->mappingFromFunction_[i];
      }

    return o;
  }


} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HXX
