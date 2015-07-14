// Copyright (C) 2015 by Félix Darricau, AIST, CNRS, EPITA
//                       Grégoire Duchemin, AIST, CNRS, EPITA
//
// This file is part of the roboptim.
//
// roboptim is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_WRAPPER_ON_MANIFOLD_HXX
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_WRAPPER_ON_MANIFOLD_HXX
# include <vector>
# include <queue>
# include <utility>
# include <iostream>
# include <typeinfo>
# include <sstream>

# include <boost/format.hpp>
# include <boost/mpl/assert.hpp>

# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>
# include <manifolds/S2.h>
namespace roboptim
{

  template <typename T>
  std::ostream&
  WrapperOnManifold<T>::print_ (std::ostream& o)
  {
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	o << (i>0?", ":"") << this->mappingFromFunction_[i];
      }

    return o;
  }

  // ---- //

  template <typename T>
  template <typename V, typename W>
  void WrapperOnManifold<T>::
  computeMapping(DescriptiveWrapper<V, W>& descWrap,
		 const mnf::Manifold& problemManifold,
		 const mnf::Manifold& functionManifold,
		 std::vector<const mnf::Manifold*> restrictedManifolds,
		 std::vector<std::pair<long, long>> restrictions)
  {
    // TODO: refactor EVERYTHING into smaller methods, and call them
    // from here. Seriously, did you check the size of this method ?

    // Assert to check the sizes of the restrictions' std::vector
    // Either they are the same, or the restrictions array is a single pair
    // If the restriction array is a single pair, it is applied to all
    // restricted manifolds.
    //
    // Also, check that the restrictions actually make sense.
    assert (restrictions.size() == 1 || (restrictedManifolds.size() == restrictions.size()));
    ROBOPTIM_DEBUG_ONLY(                                                                                    \
			for (size_t i = 0; i < restrictedManifolds.size(); ++i)                                               \
			  {                                                                                                   \
			    std::pair<long BOOST_PP_COMMA() long> restriction = restrictions[(restrictions.size() == 1?0:i)]; \
			    assert (restriction.first + restriction.second <= restrictedManifolds[i]->representationDim());   \
			  })

      // Check if each manifold appeared as much time in the restricted list as in the manifold tree
      std::map<long, int> occurenceCount;
    occurenceCount.clear();
    std::map<long, int> traversalOccurences;
    traversalOccurences.clear();
    std::map<long, const mnf::Manifold*> manifoldsById;

    // Compute the number of time each manifold should appear
    // If 0, the manifold can appear as much as it wants
    for (size_t i = 0; i < restrictedManifolds.size(); ++i)
      {
	if (occurenceCount.find(restrictedManifolds[i]->getInstanceId()) != occurenceCount.end())
	  {
	    occurenceCount[restrictedManifolds[i]->getInstanceId()] += 1;
	  }
	else
	  {
	    occurenceCount[restrictedManifolds[i]->getInstanceId()] = 1;
	  }

	if (manifoldsById.find(restrictedManifolds[i]->getInstanceId()) == manifoldsById.end())
	  {
	    manifoldsById[restrictedManifolds[i]->getInstanceId()] = restrictedManifolds[i];
	  }
      }

    // Compute the occurences of each manifold in the submanifold
    // on which the function will be instantiated
    std::function<void(const mnf::Manifold&)> checkOccurences =
      [&traversalOccurences, &checkOccurences, &manifoldsById]
      (const mnf::Manifold& manifold)
      {
	if (manifold.isElementary())
	  {
	    if (traversalOccurences.find(manifold.getInstanceId()) != traversalOccurences.end())
	      {
		traversalOccurences[manifold.getInstanceId()] += 1;
	      }
	    else
	      {
		traversalOccurences[manifold.getInstanceId()] = 1;
	      }

	    if (manifoldsById.find(manifold.getInstanceId()) == manifoldsById.end())
	      {
		manifoldsById[manifold.getInstanceId()] = &manifold;
	      }
	  }
	else
	  {
	    for (size_t i = 0; i < manifold.numberOfSubmanifolds(); ++i)
	      {
		checkOccurences(manifold(i));
	      }
	  }
      };
    checkOccurences(functionManifold);

    // Compare the two maps. If the occurenceCount value is set,
    // then the traversalOccurences value should be equal to it.

    for(std::map<long, int>::iterator ite = occurenceCount.begin(); ite != occurenceCount.end(); ++ite)
      {
	if (traversalOccurences.find(ite->first) == traversalOccurences.end())
	  {
	    std::stringstream ss;
	    ss << "Error: manifold " << manifoldsById[ite->first]->name() << " (instanceId: " << ite->first
	       << ") did not appear in the submanifold, but appeared "
	       << ite->second << " in the restriction list." << std::endl
	       << "You should not put a restriction on a manifold you are not using the function upon."
	       << std::endl;
	    throw std::runtime_error(ss.str().c_str());
	  }
	if (ite->second != traversalOccurences[ite->first])
	  {
	    std::stringstream ss;
	    ss << "Error: manifold " << manifoldsById[ite->first]->name() << " (instanceId: " << ite->first
	       << ") appeared " << traversalOccurences[ite->first] << " times in the submanifold, but only "
	       << ite->second << " in the restriction list." << std::endl
	       << "Either restrict all its occurences (even if it is restricted to the whole space) or do not restrict any."
	       << std::endl;
	    throw std::runtime_error(ss.str().c_str());
	  }
      }

    // --METHOD START-- //

    bool onTangentSpace = false;
    //	this->fct_ = descWrap.fct();

    // TODO: should be memoized for performance, although we can compute
    // an adequate map if we use a Factory pattern to create this wrapper.
    //
    // This lambda returns a std::pair of long representing a restriction
    // of the manifold. The pair contains the starting index as the first
    // element and the size of the restriction as the second element.
    std::function<std::pair<long, long>(const mnf::Manifold&, const int, int)> getRestriction =
      [&restrictedManifolds, &restrictions, &onTangentSpace]
      (const mnf::Manifold& manifold, const int searchedIndex, int currIndex)
      {
	// A (-1, -1) is equivalent to no restrictions at all
	std::pair<long, long> ans = std::make_pair(-1l, -1l);

	for (size_t i = 0; i < restrictedManifolds.size(); ++i)
	  {
	    if (manifold.getInstanceId() == restrictedManifolds[i]->getInstanceId())
	      {
		if (searchedIndex == currIndex)
		  {
		    ans = restrictions[(restrictions.size() == 1?0:i)];
		    break;
		  }
		else
		  {
		    ++currIndex;
		  }
	      }
	  }

	// If we could not find a restriction for this manifold
	// we set the restriction to the entire manifold's size
	ans.first = std::max(0l, ans.first);
	if (ans.second < 0)
	  {
	    ans.second = (onTangentSpace?manifold.tangentDim():manifold.representationDim());
	  }

	return ans;
      };

    // These helpers function will help us to track the id of the manifolds
    // we encounter while traversing the manifold trees, reusing the maps
    // we first used to check the occurences counts of manifolds
    std::function<int(const mnf::Manifold&)> getCurrentOccurenceCount =
      [&traversalOccurences]
      (const mnf::Manifold& manifold)
      {
	if (traversalOccurences.find(manifold.getInstanceId()) != traversalOccurences.end())
	  {
	    return traversalOccurences[manifold.getInstanceId()];
	  }

	return 0;
      };
    std::function<void(const mnf::Manifold&)> incrementOccurenceCount =
      [&traversalOccurences]
      (const mnf::Manifold& manifold)
      {
	if (traversalOccurences.find(manifold.getInstanceId()) != traversalOccurences.end())
	  {
	    traversalOccurences[manifold.getInstanceId()] += 1;
	  }
	else
	  {
	    traversalOccurences[manifold.getInstanceId()] = 1;
	  }
      };

    std::vector<const mnf::Manifold*> planarManifold;

    // This lambda converts a manifold tree to a std::vector of its leaf
    // which should all be elementary manifolds.
    std::function<void(const mnf::Manifold&)> manifoldTreeToVector = [&planarManifold, &manifoldTreeToVector](const mnf::Manifold& manifold)
      {
	if (manifold.isElementary())
	  {
	    planarManifold.push_back(&manifold);
	    return;
	  }

	for(size_t i = 0; i < manifold.numberOfSubmanifolds(); ++i)
	  {
	    manifoldTreeToVector(manifold(i));
	  }
      };
    manifoldTreeToVector(functionManifold);

    // This function compares a manifold tree to a std::vector
    // of elementary manifolds.
    // This is because the tree structure of the manifold is
    // not important for the mapping; only the ordering counts.
    std::function<long(const mnf::Manifold&, long)> checkManifold = [&planarManifold, &checkManifold, &getRestriction, &getCurrentOccurenceCount, &incrementOccurenceCount](const mnf::Manifold& manifold, long index)
      {
	if (manifold.isElementary())
	  {
	    const mnf::Manifold* myManifold = planarManifold[static_cast<size_t>(index)];
	    bool sameType = myManifold->getTypeId() == manifold.getTypeId();
	    bool isNotRealSpace = myManifold->getTypeId() != mnf::RealSpace(1).getTypeId();
	    long size1 = getRestriction(*myManifold, getCurrentOccurenceCount(*myManifold), 0).second;
	    long size2 = manifold.representationDim();
	    bool sameSize = size1 == size2;
	    incrementOccurenceCount(*myManifold);

	    if (sameType && (isNotRealSpace || sameSize))
	      {
		return ++index;
	      }

	    if (sameType && !isNotRealSpace && !sameSize)
	      {
		std::cerr << "RealSpaces must be of same size: Expecting "
		<< size1 << " but got " << size2 << " instead." << std::endl;
	      }
	    return -1l;
	  }

	for (long i = 0; i < static_cast<long>(manifold.numberOfSubmanifolds()); ++i)
	  {
	    index = checkManifold(manifold(static_cast<size_t>(i)), index);
	    if (index < 0)
	      {
		std::cerr << "A child has failed" << std::endl;
		return -1l;
	      }
	  }

	return index;
      };

    // Clear the traversalOccurences map to ensure a clean use of
    // the getCurrentOccurenceCount and incrementOccurenceCount
    // helper lambda functions.
    traversalOccurences.clear();

    // Raise an exception if the numbers of matched manifolds is not
    // equal to the total number of manifolds to match.
    long manifoldsMatched = checkManifold(descWrap.manifold(), 0);

    if (manifoldsMatched != static_cast<long>(planarManifold.size()))
      {
	throw std::runtime_error("ERROR: instantiation and descriptive manifolds are not equivalent");
      }

    // This lambda computes the dimension of the input space while taking the restrictions into account
    std::function<long(const mnf::Manifold&)> computeRestrictedDimension = [&computeRestrictedDimension, &getRestriction, &getCurrentOccurenceCount, &incrementOccurenceCount](const mnf::Manifold& manifold)
      {
	long mySize = 0;

	if (manifold.isElementary())
	  {
	    mySize = getRestriction(manifold, getCurrentOccurenceCount(manifold), 0).second;
	    incrementOccurenceCount(manifold);
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

    // Clear the traversalOccurences map to ensure a clean use of
    // the getCurrentOccurenceCount and incrementOccurenceCount
    // helper lambda functions.
    traversalOccurences.clear();
    this->mappingFromFunctionSize_ = computeRestrictedDimension(functionManifold);
    this->mappingFromFunction_ = new size_t[this->mappingFromFunctionSize_];

    // Clear the traversalOccurences map to ensure a clean use of
    // the getCurrentOccurenceCount and incrementOccurenceCount
    // helper lambda functions.
    traversalOccurences.clear();
    onTangentSpace = true;
    this->tangentMappingFromFunctionSize_ = computeRestrictedDimension(functionManifold);
    this->tangentMappingFromFunction_ = new size_t[this->tangentMappingFromFunctionSize_];
    onTangentSpace = false;

    this->mappedInput_ = Eigen::VectorXd::Zero(this->mappingFromFunctionSize_);

    // This lambda computes the actual mapping between a manifold and the one
    // in its place in the global manifold of the problem.
    std::function<long(const mnf::Manifold&, long, int, long, long)> getStartingIndexOfManifold = [&getStartingIndexOfManifold, this, &getRestriction, &onTangentSpace](const mnf::Manifold& manifold, long targetId, int targetIndex, long functionStartIndex, long startIndex)
      {
	if (targetId == manifold.getInstanceId())
	  {
	    // We found the manifold
	    // Write its indexes and exit
	    std::pair<long, long> restriction = getRestriction(manifold, targetIndex, 0);

	    for (long i = 0; i < restriction.second; ++i)
	      {
		if (onTangentSpace)
		  {
		    this->tangentMappingFromFunction_[static_cast<size_t>(functionStartIndex + i)]
		      = static_cast<size_t>(startIndex + i + restriction.first);
		  }
		else
		  {
		    this->mappingFromFunction_[static_cast<size_t>(functionStartIndex + i)]
		      = static_cast<size_t>(startIndex + i + restriction.first);
		  }
	      }

	    return -1l;
	  }
	else
	  {
	    if (manifold.isElementary())
	      {
		return startIndex + (onTangentSpace?manifold.tangentDim():manifold.representationDim());
	      }
	    else
	      {
		for (size_t i = 0; i < manifold.numberOfSubmanifolds() && startIndex >= 0; ++i)
		  {
		    startIndex = getStartingIndexOfManifold(manifold(i), targetId, targetIndex, functionStartIndex, startIndex);
		  }

		return startIndex;
	      }
	  }
      };

    // This lambda recursively traverse the manifold tree, computing the mapping
    /// for each leaf (elementary manifold) it reaches.
    std::function<int(const mnf::Manifold&, int)> traverseFunctionManifold = [&traverseFunctionManifold, &getStartingIndexOfManifold, &problemManifold, &getRestriction, &onTangentSpace, &getCurrentOccurenceCount, &incrementOccurenceCount](const mnf::Manifold& manifold, int startIndex)
      {
	if (manifold.isElementary())
	  {
	    int targetIndex = getCurrentOccurenceCount(manifold);
	    incrementOccurenceCount(manifold);
	    getStartingIndexOfManifold(problemManifold, manifold.getInstanceId(), targetIndex, startIndex, 0);
	    int retValue = static_cast<int>(startIndex + getRestriction(manifold, targetIndex, 0).second);

	    return retValue;
	  }

	for (size_t i = 0; i < manifold.numberOfSubmanifolds() && startIndex >= 0; ++i)
	  {
	    startIndex = traverseFunctionManifold(manifold(i), startIndex);
	  }

	return static_cast<int>(startIndex);
      };

    // Clear the traversalOccurences map to ensure a clean use of
    // the getCurrentOccurenceCount and incrementOccurenceCount
    // helper lambda functions.
    traversalOccurences.clear();
    // Computes the mapping
    traverseFunctionManifold(functionManifold, 0);

    // Clear the traversalOccurences map to ensure a clean use of
    // the getCurrentOccurenceCount and incrementOccurenceCount
    // helper lambda functions.
    traversalOccurences.clear();
    onTangentSpace = true;
    traverseFunctionManifold(functionManifold, 0);
  }

  template <typename T>
  void
  WrapperOnManifold<T>::mapArgument (const_argument_ref argument)
    const
  {
    for (long i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	this->mappedInput_(i) = argument(static_cast<long>(this->mappingFromFunction_[i]));
      }
  }

  template <typename T>
  void
  WrapperOnManifold<T>::impl_compute
  (result_ref result, const_argument_ref x)
    const
  {
    this->mapArgument(x);
    (*this->fct_)(result, this->mappedInput_);
  }

  template <typename T>
  WrapperOnManifold<T>::~WrapperOnManifold()

  {
    delete [] this->mappingFromFunction_;
    delete [] this->tangentMappingFromFunction_;
  }

  template <typename T>
  void
  WrapperOnManifold<T>::unmapTangentJacobian(mnf::RefMat jacobian)
    const
  {
    for (long i = 0; i < this->tangentMappingFromFunctionSize_; ++i)
      {
	jacobian.col(static_cast<long>(this->tangentMappingFromFunction_[i])) = this->tangentMappedJacobian.col(i);
      }
  }

  template <typename T>
  void
  WrapperOnManifold<T>::unmapGradient(gradient_ref gradient) const
  {
    assert(gradient.cols() == this->inputSize());

    for (int i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	gradient.coeffRef(static_cast<size_type>(this->mappingFromFunction_[i])) = this->mappedGradient_.coeffRef(i);
      }
  }

  template <typename T>
  void
  WrapperOnManifold<T>::impl_gradient (gradient_ref gradient,
					const_argument_ref argument,
					size_type functionId)
    const
  {
    this->mapArgument(argument);
    this->fctDiff_->gradient(this->mappedGradient_, this->mappedInput_, functionId);

    gradient.setZero();
    unmapGradient(gradient);
  }

  template <typename T>
  void
  WrapperOnManifold<T>::impl_jacobian (jacobian_ref jacobian,
					const_argument_ref argument)
    const
  {
    this->mapArgument(argument);
    jacobian.setZero();

    this->fctDiff_->jacobian(this->mappedJacobian_, this->mappedInput_);
    this->unmapJacobian(jacobian);
  }

  template <typename T>
  void
  WrapperOnManifold<T>::manifold_jacobian (mnf::RefMat jacobian,
					    const_argument_ref argument)
    const
  {
    this->mapArgument(argument);
    jacobian.setZero();
    this->mappedJacobian_.setZero();
    this->tangentMappedJacobian.setZero();

    this->fctDiff_->jacobian(this->mappedJacobian_, this->mappedInput_);

    this->applyDiff();

    this->unmapTangentJacobian(jacobian);
  }

  template<typename T>
  void
  WrapperOnManifold<T>::unmapHessian(hessian_ref hessian) const
  {
    assert (hessian.rows() == this->inputSize());
    assert (hessian.cols() == this->inputSize());

    for (int i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	for (int j = 0; j < this->inputSize(); ++j)
	  hessian.coeffRef(i,j) = this->mappedHessian_.coeffRef(i,j);
      }
  }

  template<typename T>
  void
  WrapperOnManifold<T>::impl_hessian(hessian_ref hessian, const_argument_ref argument, size_type functionId) const
  {
    this->mapArgument(argument);
    hessian.setZero();

    this->fctTwiceDiff_->hessian(this->mappedHessian_, this->mappedInput_, functionId);
    unmapHessian(hessian);
  }

  template<typename T>
  void
  WrapperOnManifold<T>::unmapJacobian(jacobian_ref jacobian) const
  {
    assert(jacobian.cols() == this->inputSize());
    assert(jacobian.rows() == this->outputSize());

    for (int i = 0; i < this->mappingFromFunctionSize_; ++i)
      {
	jacobian.col(static_cast<int>(this->mappingFromFunction_[i])) = this->mappedJacobian_.col(i);
      }
  }

  template<>
  inline
  void
  WrapperOnManifold<EigenMatrixDense>::applyDiff() const
  {
    this->manifold_->applyDiffRetractation(this->tangentMappedJacobian, this->mappedJacobian_, this->mappedInput_);
  }

  template<>
  inline
  void
  WrapperOnManifold<EigenMatrixSparse>::applyDiff() const
  {
    if (this->manifold_->getTypeId() != mnf::RealSpace(1).getTypeId())
      throw std::runtime_error("No Sparse support for non RealSpaces manifolds");
    else
      this->tangentMappedJacobian = this->mappedJacobian_;
  }

  template <typename T>
  const mnf::Manifold*
  WrapperOnManifold<T>::getManifold() const
  {
    return this->manifold_;
  }

  template <typename T>
  typename WrapperOnManifold<T>::flag_t
  WrapperOnManifold<T>::getFlags() const
  {
    return flags|fct_->getFlags();
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_WRAPPER_ON_MANIFOLD_HXX
