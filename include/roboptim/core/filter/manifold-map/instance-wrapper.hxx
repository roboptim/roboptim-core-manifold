#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
# include <boost/format.hpp>
# include <vector>
# include <queue>
# include <manifolds/Manifold.h>
# include <iostream>

namespace roboptim
{

  template <typename U>
  InstanceWrapper<U>::InstanceWrapper
  (boost::shared_ptr<U> origin,
   pgs::Manifold& problemManifold,
   pgs::Manifold& functionManifold)
    : detail::AutopromoteTrait<U>::T_type
      (origin->inputSize (),
       origin->outputSize (),
       (boost::format ("%1%")
	% origin->getName ()).str ()),
      origin_ (origin)
  {
    this->mappingFromProblemSize_ = problemManifold.representationDim();
    this->mappingFromProblem_ = new int[this->mappingFromProblemSize_];
    this->mappingFromFunctionSize_ = functionManifold.representationDim();
    this->mappingFromFunction_ = new int[this->mappingFromFunctionSize_];

    std::cout << "this->mappingFromProblemSize_: " << this->mappingFromProblemSize_ << std::endl;
    std::cout << "this->mappingFromFunctionSize_: " << this->mappingFromFunctionSize_ << std::endl;

    std::function<int(const pgs::Manifold&, std::string, int, int)> getStartingIndexOfManifold = [&getStartingIndexOfManifold, this](const pgs::Manifold& manifold, std::string targetName, int functionStartIndex, int startIndex)
      {
	if (!targetName.compare(manifold.name()))
	  {
	    // We found the manifold
	    // Write its index and exit
	    std::cout << "Found a match for " << manifold.name()  << std::endl;
	    std::cout << "startIndex: " << startIndex << std::endl;
	    std::cout << "functionStartIndex: " << functionStartIndex << std::endl;

	    for (size_t i = 0; i < manifold.representationDim(); ++i)
	      {
		std::cout << "i: " << i << std::endl;
		this->mappingFromProblem_[startIndex + i] = static_cast<int>(functionStartIndex + i);
		this->mappingFromFunction_[startIndex + i] = static_cast<int>(startIndex + i);
	      }

	    return -1;
	  }
	else
	  {
	    if (manifold.isElementary())
	      {
		return static_cast<int>(startIndex + manifold.representationDim());
	      }
	    else
	      {
		for (size_t i = 0; i < manifold.numberOfSubmanifolds() && startIndex >= 0; ++i)
		  {
		    startIndex = getStartingIndexOfManifold(manifold(i), targetName, functionStartIndex, startIndex);
		  }

		return static_cast<int>(startIndex);
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
    delete this->mappingFromProblem_;
    delete this->mappingFromFunction_;
  }

  // FIXME: temp, complete those methods
  template <typename U>
  void
  InstanceWrapper<U>::impl_compute
  (result_ref result, const_argument_ref x)
    const
  {
    std::cout << "result: " << result << std::endl;
    std::cout << "x: " << x << std::endl;
  }

  template <typename U>
  void
  InstanceWrapper<U>::impl_gradient (gradient_ref gradient,
			 const_argument_ref argument,
			 size_type functionId)
    const
  {
    std::cout << "gradient: " << std::endl << gradient << std::endl;
    std::cout << "argument: " << std::endl << argument << std::endl;
    std::cout << "functionId: " << std::endl << functionId << std::endl;
  }

  template <typename U>
  void
  InstanceWrapper<U>::impl_jacobian (jacobian_ref jacobian,
			 const_argument_ref argument)
    const
  {
    std::cout << "jacobian: " << std::endl << jacobian << std::endl;
    std::cout << "argument: " << std::endl << argument << std::endl;
  }

  template <typename U>
  std::ostream&
  InstanceWrapper<U>::print (std::ostream& o)
  {
    for (size_t i = 0; i < this->mappingFromProblemSize_; ++i)
      {
	o << (i>0?", ":"") << this->mappingFromProblem_[i];
      }
    o << "\n";
    for (size_t i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	o << (i>0?", ":"") << this->mappingFromFunction_[i];
      }

    return o;
  }


} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
