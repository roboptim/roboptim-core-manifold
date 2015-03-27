#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
# include <vector>
# include <queue>
# include <utility>
# include <iostream>

# include <boost/format.hpp>
# include <boost/mpl/assert.hpp>

# include <roboptim/core/filter/manifold-map/descriptive-wrapper.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>
# include <manifolds/S2.h>
namespace roboptim
{

  template <typename U>
  InstanceWrapper<U>::InstanceWrapper
  (boost::shared_ptr<DescriptiveWrapper<U>> descWrap,
   const pgs::Manifold& problemManifold,
   const pgs::Manifold& functionManifold,
   std::vector<const pgs::Manifold*> restrictedManifolds,
   std::vector<std::pair<long, long>> restrictions)
    : detail::AutopromoteTrait<U>::T_type
      (problemManifold.representationDim(),
       descWrap->fct().outputSize (),
       (boost::format ("%1%")
	% descWrap->fct().getName ()).str ()),
      descWrap_ (descWrap)
  {
    // Assert to check the sizes of the restrictions' std::vector
    // Either they are the same, or the restrictions array is a single pair
    // If the restriction array is a single pair, it is applied to all
    // restricted manifolds.
    //
    // Also, check that the restrictions actually make sense.
    assert (restrictions.size() == 1 || (restrictedManifolds.size() == restrictions.size()));
    ROBOPTIM_DEBUG_ONLY(\
      for (size_t i = 0; i < restrictedManifolds.size(); ++i)\
	{\
	  std::pair<long BOOST_PP_COMMA() long> restriction = restrictions[(restrictions.size() == 1?0:i)];\
	  assert (restriction.first + restriction.second <= restrictedManifolds[i]->representationDim());\
	})


    // TODO: should be memoized for performance, although we can compute
    // an adequate map if we use a Factory pattern to create this wrapper.
    std::function<std::pair<long, long>(const pgs::Manifold&)> getRestriction = [&restrictedManifolds, &restrictions, &getRestriction](const pgs::Manifold& manifold)
      {
	// A (-1, -1) is equivalent to no restrictions at all
	std::pair<long, long> ans = std::make_pair(-1l, -1l);

	for (size_t i = 0; i < restrictedManifolds.size(); ++i)
	  {
	    // TODO: replace by an id check
	    if (!manifold.name().compare(restrictedManifolds[i]->name()))
	      {
		ans = restrictions[(restrictions.size() == 1?0:i)];
		break;
	      }
	  }

	ans.first = std::max(0l, ans.first);
	if (ans.second < 0)
	  {
	    ans.second = manifold.representationDim();
	  }

	return ans;
      };

      std::function<bool(const pgs::Manifold&, const pgs::Manifold&)> checkManifoldMatch = [&checkManifoldMatch, &getRestriction](const pgs::Manifold& fManifold, const pgs::Manifold& dManifold)
      {
	if (dManifold.isElementary() != fManifold.isElementary())
	  {
	    // Both manifold should be either elementary or composed at this point
	    std::cout << "NONSENSE" << std::endl;
	    std::cout << "fManifold.name(): " << fManifold.name() << std::endl;
	    return false;
	  }

	if (dManifold.isElementary())
	  {
	    bool sameType = !std::strcmp(typeid(fManifold).name(), typeid(dManifold).name());
	    bool isNotRealSpace = std::strcmp(typeid(fManifold).name(), typeid(pgs::RealSpace).name());
	    bool sameSize = getRestriction(fManifold).second == dManifold.representationDim();


	    return sameType && (isNotRealSpace || sameSize);
	  }

	// Both manifolds are composed
	if (dManifold.numberOfSubmanifolds() != fManifold.numberOfSubmanifolds())
	  {
	    std::cout << "Not enough children" << std::endl;
	    return false;
	  }

	for (size_t i = 0; i < dManifold.numberOfSubmanifolds(); ++i)
	  {
	    if (!checkManifoldMatch(fManifold(i), dManifold(i)))
	      {
		std::cout << "Children do not match" << std::endl;
		return false;
	      }
	  }

	return true;
      };

      if (!checkManifoldMatch(functionManifold, descWrap_->manifold()))
	{
	  throw std::runtime_error("LOL");
	}

    std::function<long(const pgs::Manifold&)> computeRestrictedDimension = [&computeRestrictedDimension, &getRestriction](const pgs::Manifold& manifold)
      {
	long mySize = 0;

	if (manifold.isElementary())
	  {
	    mySize = getRestriction(manifold).second;
	  }
	else
	  {
	    for (size_t i = 0; i < manifold.numberOfSubmanifolds(); ++i)
	      {
		mySize += computeRestrictedDimension(manifold(i));
	      }
	  }

	return mySize;
      };

    this->mappingFromFunctionSize_ = computeRestrictedDimension(functionManifold);
    this->mappingFromFunction_ = new size_t[this->mappingFromFunctionSize_];

    this->mappedInput_ = Eigen::VectorXd::Zero(this->mappingFromFunctionSize_);
    this->mappedGradient_ = Eigen::VectorXd::Zero(this->mappingFromFunctionSize_);
    this->mappedJacobian_ = Eigen::MatrixXd::Zero(descWrap->fct().outputSize(), this->mappingFromFunctionSize_);


    std::function<long(const pgs::Manifold&, std::string, long, long)> getStartingIndexOfManifold = [&getStartingIndexOfManifold, this, &getRestriction](const pgs::Manifold& manifold, std::string targetName, long functionStartIndex, long startIndex)
      {
	// TODO: replace by an id check
	if (!targetName.compare(manifold.name()))
	  {
	    // We found the manifold
	    // Write its indexes and exit
	    std::pair<long, long> restriction = getRestriction(manifold);

	    for (long i = 0; i < restriction.second; ++i)
	      {
		this->mappingFromFunction_[static_cast<size_t>(functionStartIndex + i)]
		  = static_cast<size_t>(startIndex + i + restriction.first);
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

    std::function<int(const pgs::Manifold&, int)> traverseFunctionManifold = [&traverseFunctionManifold, &getStartingIndexOfManifold, &problemManifold, &getRestriction](const pgs::Manifold& manifold, int startIndex)
      {
	if (manifold.isElementary())
	  {
	    getStartingIndexOfManifold(problemManifold, manifold.name(), startIndex, 0);
	    return static_cast<int>(startIndex + getRestriction(manifold).second);
	  }

	for (size_t i = 0; i < manifold.numberOfSubmanifolds() && startIndex >= 0; ++i)
	  {
	    startIndex = traverseFunctionManifold(manifold(i), startIndex);
	  }

	return static_cast<int>(startIndex);
      };

    // Computes the mapping
    traverseFunctionManifold(functionManifold, 0);

  }

  template <typename U>
  InstanceWrapper<U>::~InstanceWrapper()
  {
    delete [] this->mappingFromFunction_;
  }

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
  InstanceWrapper<U>::unmapGradient(gradient_ref gradient)
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

    for (long j = 0; j < jacobian.rows(); ++j)
      {
	descWrap_->fct().gradient(this->mappedGradient_, this->mappedInput_, j);
	this->unmapGradient(jacobian.row(j));
      }

  }

  template <typename U>
  std::ostream&
  InstanceWrapper<U>::print_ (std::ostream& o)
  {
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	o << (i>0?", ":"") << this->mappingFromFunction_[i];
      }

    return o;
  }


} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_INSTANCE_WRAPPER_HXX
