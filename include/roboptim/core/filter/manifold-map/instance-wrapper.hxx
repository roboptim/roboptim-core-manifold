#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
# include <boost/format.hpp>
# include <vector>
# include <queue>
# include <manifolds/Manifold.h>
# include <iostream>

# include <roboptim/core/filter/manifold-map/descriptive-wrapper.hh>

namespace roboptim
{

  template <typename U>
  InstanceWrapper<U>::InstanceWrapper
  (boost::shared_ptr<DescriptiveWrapper<U>> descWrap,
     pgs::Manifold& problemManifold,
     pgs::Manifold& functionManifold)
    : detail::AutopromoteTrait<U>::T_type
      (descWrap->fct().inputSize (),
       descWrap->fct().outputSize (),
       (boost::format ("%1%")
	% descWrap->fct().getName ()).str ()),
      descWrap_ (descWrap)
  {
    this->mappingFromProblemSize_ = problemManifold.representationDim();
    this->mappingFromProblem_ = new size_t[this->mappingFromProblemSize_];
    this->mappingFromFunctionSize_ = functionManifold.representationDim();
    this->mappingFromFunction_ = new size_t[this->mappingFromFunctionSize_];

    std::cout << "this->mappingFromProblemSize_: " << this->mappingFromProblemSize_ << std::endl;
    std::cout << "this->mappingFromFunctionSize_: " << this->mappingFromFunctionSize_ << std::endl;

    std::function<long(const pgs::Manifold&, std::string, long, long)> getStartingIndexOfManifold = [&getStartingIndexOfManifold, this](const pgs::Manifold& manifold, std::string targetName, long functionStartIndex, long startIndex)
      {
	if (!targetName.compare(manifold.name()))
	  {
	    // We found the manifold
	    // Write its index and exit
	    std::cout << "Found a match for " << manifold.name()  << std::endl;
	    std::cout << "startIndex: " << startIndex << std::endl;
	    std::cout << "functionStartIndex: " << functionStartIndex << std::endl;

	    for (long i = 0; i < manifold.representationDim(); ++i)
	      {
		std::cout << "i: " << i << std::endl;
		this->mappingFromProblem_[static_cast<size_t>(startIndex + i)] = static_cast<size_t>(functionStartIndex + i);
		this->mappingFromFunction_[static_cast<size_t>(functionStartIndex + i)] = static_cast<size_t>(startIndex + i);
	      }

	    return -1l;
	  }
	else
	  {
	    if (manifold.isElementary())
	      {
		return startIndex + manifold.representationDim();
	      }
	    else
	      {
		for (size_t i = 0; i < manifold.numberOfSubmanifolds() && startIndex >= 0; ++i)
		  {
		    startIndex = getStartingIndexOfManifold(manifold(i), targetName, functionStartIndex, startIndex);
		  }

		return startIndex;
	      }
	  }
      };

    std::function<int(const pgs::Manifold&, int)> traverseFunctionManifold = [&traverseFunctionManifold, &getStartingIndexOfManifold, &problemManifold](const pgs::Manifold& manifold, int startIndex)
      {
	if (manifold.isElementary())
	  {
	    std::cout << "Matching function manifold: " << manifold.name()  << std::endl;
	    getStartingIndexOfManifold(problemManifold, manifold.name(), startIndex, 0);
	    return static_cast<int>(startIndex + manifold.representationDim());
	  }

	for (size_t i = 0; i < manifold.numberOfSubmanifolds() && startIndex >= 0; ++i)
	  {
	    startIndex = traverseFunctionManifold(manifold(i), startIndex);
	  }

	return static_cast<int>(startIndex);
      };

    traverseFunctionManifold(functionManifold, 0);
  }

  template <typename U>
  InstanceWrapper<U>::~InstanceWrapper()
  {
    delete [] this->mappingFromProblem_;
    delete [] this->mappingFromFunction_;
  }

  // FIXME: temp, complete those methods
  template <typename U>
  void
  InstanceWrapper<U>::mapArgument (const_argument_ref argument)
    const
  {
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	this->mappedInput_(i) = argument(static_cast<long>(this->mappingFromFunction_[i]));
      }
  }

  template <typename U>
  void
  InstanceWrapper<U>::unmapGradient(gradient_t gradient)
    const
  {
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	gradient(static_cast<long>(this->mappingFromFunction_[i])) = this->mappedGradient_(i);
      }
  }

  template <typename U>
  void
  InstanceWrapper<U>::impl_compute
  (result_ref result, const_argument_ref x)
    const
  {
    this->mapArgument(x);

    descWrap_->fct()(result, this->mappedInput_);
  }

  template <typename U>
  void
  InstanceWrapper<U>::impl_gradient (gradient_ref gradient,
			 const_argument_ref argument,
			 size_type functionId)
    const
  {
    this->mapArgument(argument);
    descWrap_->fct().gradient(this->mappedGradient_, this->mappedInput_, functionId);

    gradient.setZero();
    this->unmapGradient(gradient);
  }

  template <typename U>
  void
  InstanceWrapper<U>::impl_jacobian (jacobian_ref jacobian,
			 const_argument_ref argument)
    const
  {
    this->mapArgument(argument);
    jacobian.setZero();

    for (long j = 0; j < jacobian.cols(); ++j)
      {
	descWrap_->fct().gradient(this->mappedGradient_, this->mappedInput_, j);
	this->unmapGradient(jacobian.col(j));
      }

  }

  template <typename U>
  std::ostream&
  InstanceWrapper<U>::print_ (std::ostream& o)
  {
    for (long i = 0; i < this->mappingFromProblemSize_; ++i)
      {
	o << (i>0?", ":"") << this->mappingFromProblem_[i];
      }
    o << "\n";
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	o << (i>0?", ":"") << this->mappingFromFunction_[i];
      }
    o << "\n";

    return o;
  }


} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
