// Copyright (C) 2015 by Grégoire Duchemin, AIST, CNRS, EPITA
//                       Félix Darricau, AIST, CNRS, EPITA
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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
# include <vector>
# include <iostream>
# include <utility>
# include <boost/shared_ptr.hpp>
# include <boost/make_shared.hpp>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>
# include <manifolds/S2.h>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  /// \brief Maps a DescriptiveWrapper to a instance of a submanifold.
  ///
  /// \tparam U input roboptim function type.
  template <typename U>
  class FunctionOnManifold : public detail::AutopromoteTrait<U>::T_type
  {
  public:
    typedef typename detail::AutopromoteTrait<U>::T_type parentType_t;
    ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (parentType_t);

    typedef boost::shared_ptr<FunctionOnManifold> FunctionOnManifoldShPtr_t;

    /// \brief Creates the mapping between the manifolds.
    ///
    /// Programmers should ideally not directly use this constructor, amd
    /// rely on the macros defined in the file manifold-map.hh to get the well
    /// defined FunctionOnManifold type, ready to be instantiated.
    ///
    /// \tparam V input function type templating the DescriptiveWrapper.
    /// \tparam W descriptive manifold type templating the DescriptiveWrapper.
    ///
    /// \param descwrap instance of the DescriptiveWrapper
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    /// \param restrictedManifolds a list of elementary Manifolds to be restricted to a part of themselves
    /// \param restrictions the restrictions applying to the selected manifolds, represented as (startingIndex, size). If a single one is given, it will apply to all restricted manifolds.
    template <typename V, typename W>
    explicit FunctionOnManifold
    (boost::shared_ptr<DescriptiveWrapper<V, W>> descWrap,
     const pgs::Manifold& problemManifold,
     const pgs::Manifold& functionManifold,
     std::vector<const pgs::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
    : detail::AutopromoteTrait<U>::T_type
    (problemManifold.representationDim(),
     descWrap->fct().outputSize (),
     (boost::format ("%1%")
      % descWrap->fct().getName ()).str ()),
      fct_(descWrap->fctPointer()),
      manifold_(descWrap->manifoldPointer())
    {
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

      bool onTangentSpace = false;
	//	this->fct_ = descWrap->fct();

	// TODO: should be memoized for performance, although we can compute
	// an adequate map if we use a Factory pattern to create this wrapper.
	//
	// This lambda returns a std::pair of long representing a restriction
	// of the manifold. The pair contains the starting index as the first
	// element and the size of the restriction as the second element.
	std::function<std::pair<long, long>(const pgs::Manifold&)> getRestriction = [&restrictedManifolds, &restrictions, &getRestriction, &onTangentSpace](const pgs::Manifold& manifold)
	{
	  // A (-1, -1) is equivalent to no restrictions at all
	  std::pair<long, long> ans = std::make_pair(-1l, -1l);

	  for (size_t i = 0; i < restrictedManifolds.size(); ++i)
	    {
	      if (manifold.getInstanceId() == restrictedManifolds[i]->getInstanceId())
		{
		  ans = restrictions[(restrictions.size() == 1?0:i)];
		  break;
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

      std::vector<const pgs::Manifold*> planarManifold;

      // This lambda converts a manifold tree to a std::vector of its leaf
      // which should all be elementary manifolds.
      std::function<void(const pgs::Manifold&)> manifoldTreeToVector = [&planarManifold, &manifoldTreeToVector](const pgs::Manifold& manifold)
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
      std::function<long(const pgs::Manifold&, long)> checkManifold = [&planarManifold, &checkManifold, &getRestriction](const pgs::Manifold& manifold, long index)
	{
	  if (manifold.isElementary())
	    {
	      bool sameType = planarManifold[static_cast<size_t>(index)]->getTypeId() == manifold.getTypeId();
	      bool isNotRealSpace = planarManifold[static_cast<size_t>(index)]->getTypeId() != pgs::RealSpace(1).getTypeId();
	      bool sameSize = getRestriction((*planarManifold[static_cast<size_t>(index)])).second == manifold.representationDim();

	      if (sameType && (isNotRealSpace || sameSize))
		{
		  return ++index;
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

      // Raise an exception if the numbers of matched manifolds is not
      // equal to the total number of manifolds to match.
      long manifoldsMatched = checkManifold(descWrap->manifold(), 0);

      if (manifoldsMatched != static_cast<long>(planarManifold.size()))
	{
	  throw std::runtime_error("ERROR: instantiation and descriptive manifolds are not equivalent");
	}

      // This lambda computes the dimension of the input space while taking the restrictions into account
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

      onTangentSpace = true;
      this->tangentMappingFromFunctionSize_ = computeRestrictedDimension(functionManifold);
      this->tangentMappingFromFunction_ = new size_t[this->tangentMappingFromFunctionSize_];
      onTangentSpace = false;

      this->mappedInput_ = Eigen::VectorXd::Zero(this->mappingFromFunctionSize_);
      this->mappedGradient_ = Eigen::VectorXd::Zero(this->mappingFromFunctionSize_);
      this->mappedJacobian_ = Eigen::MatrixXd::Zero(descWrap->fct().outputSize(), this->mappingFromFunctionSize_);
      this->tangentMappedJacobian = Eigen::MatrixXd::Zero(descWrap->fct().outputSize(), this->tangentMappingFromFunctionSize_);

      // This lambda computes the actual mapping between a manifold and the one
      // in its place in the global manifold of the problem.
      std::function<long(const pgs::Manifold&, long, long, long)> getStartingIndexOfManifold = [&getStartingIndexOfManifold, this, &getRestriction, &onTangentSpace](const pgs::Manifold& manifold, long targetId, long functionStartIndex, long startIndex)
	{
	  if (targetId == manifold.getInstanceId())
	    {
	      // We found the manifold
	      // Write its indexes and exit
	      std::pair<long, long> restriction = getRestriction(manifold);

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
		      startIndex = getStartingIndexOfManifold(manifold(i), targetId, functionStartIndex, startIndex);
		    }

		  return startIndex;
		}
	    }
	};

      // This lambda recursively traverse the manifold tree, computing the mapping
      // for each leaf (elementary manifold) it reaches.
      std::function<int(const pgs::Manifold&, int)> traverseFunctionManifold = [&traverseFunctionManifold, &getStartingIndexOfManifold, &problemManifold, &getRestriction, &onTangentSpace](const pgs::Manifold& manifold, int startIndex)
	{
	  if (manifold.isElementary())
	    {
	      getStartingIndexOfManifold(problemManifold, manifold.getInstanceId(), startIndex, 0);
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

      onTangentSpace = true;
      traverseFunctionManifold(functionManifold, 0);
    }


    /// \brief Creates the mapping between the manifolds without restrictions.
    ///
    /// Programmers should ideally not directly use this constructor, amd
    /// rely on the macros defined in the file manifold-map.hh to get the well
    /// defined FunctionOnManifold type, ready to be instantiated.
    ///
    /// \tparam V input function type templating the DescriptiveWrapper.
    /// \tparam W descriptive manifold type templating the DescriptiveWrapper.
    ///
    /// \param fct instance of the DescriptiveWrapper
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    template<typename V, typename W>
    explicit FunctionOnManifold (boost::shared_ptr<DescriptiveWrapper<V, W>> fct,
			      const pgs::Manifold& problemManifold,
			      const pgs::Manifold& functionManifold)
      : FunctionOnManifold(fct, problemManifold, functionManifold,
	     std::vector<const pgs::Manifold*>(), std::vector<std::pair<long, long>>())
    {
    }

    ~FunctionOnManifold ();

    void impl_compute (result_ref result, const_argument_ref x)
      const;

    void impl_gradient (gradient_ref gradient,
			const_argument_ref argument,
			size_type functionId = 0)
      const;
    void impl_jacobian (jacobian_ref jacobian,
			const_argument_ref arg)
      const;

    void manifold_jacobian (pgs::RefMat jacobian,
			    const_argument_ref arg)
      const;

    std::ostream& print_(std::ostream& o);
  private:
  public:
    /// \brief the function.
    boost::shared_ptr<U> fct_;
    /// \brief the problem manifold.
    boost::shared_ptr<pgs::Manifold> manifold_;

    /// \brief array representing the restricted mapping
    size_t* mappingFromFunction_;
    /// \brief size of the array
    long mappingFromFunctionSize_;

    /// \brief array representing the restricted mapping for the tangent problem
    size_t* tangentMappingFromFunction_;
    /// \brief size of the array
    long tangentMappingFromFunctionSize_;

    /// \brief new input mapped to the restricted problem
    mutable vector_t mappedInput_;
    /// \brief new gradient mapped to the restricted problem
    mutable gradient_t mappedGradient_;
    /// \brief new jacobian mapped to the restricted problem
    mutable jacobian_t mappedJacobian_;

    // Dirty copy ONLY
    mutable Eigen::MatrixXd tangentMappedJacobian;

    /// \brief map the input to the restricted problem
    ///
    /// \param argument the argument to map
    void mapArgument(const_argument_ref argument)
      const;

    /// \brief gets the gradient from the restricted problem
    ///
    /// \param gradient the output gradient
    /// \param mappedGradient the gradient computed to unmap
    void unmapGradient(gradient_ref gradient, const_gradient_ref mappedGradient)
      const;

    /// \brief unmap the jacobian from the restricted problem
    ///
    /// \param jacobian the jacobian to unmap
    void unmapTangentJacobian(pgs::RefMat jacobian)
      const;
  };

  template <typename U>
  std::ostream&
  operator<<(std::ostream& o, FunctionOnManifold<U>& instWrap)
  {
    return instWrap.print_(o);
  }

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/manifold-map/decorator/function-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
