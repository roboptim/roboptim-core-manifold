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
#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_PROBLEM_FACTORY_HXX
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_PROBLEM_FACTORY_HXX

# include <memory>

# include <boost/pointer_cast.hpp>

# include <roboptim/core/debug.hh>
# include <roboptim/core/manifold-map/util.hh>

namespace roboptim {

template<typename T>
void ManifoldProblemFactory<T>::addElementaryManifolds(const mnf::Manifold& instanceManifold)
{
  constraintsManifold_.addManifold(instanceManifold);
}

// ---- //

template<typename T>
template<class V, class W>
BoundsAndScalingSetter<T> ManifoldProblemFactory<T>::addConstraint(std::shared_ptr<DescriptiveWrapper<V, W>> descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions)
{
  this->addElementaryManifolds(instanceManifold);
  this->boundsAndScaling_.push_back(std::make_pair(typename V::intervals_t(), typename Problem<T>::scaling_t()));

  std::pair<typename GenericFunction<T>::intervals_t, typename Problem<T>::scaling_t>* bNSPair = &(this->boundsAndScaling_.back());

  // We need the type of the function to instantiate the WrapperOnManifold
  // wrapper, hence why we capture the information we need here and wait
  // for the global manifold to be defined to execute this lambda.
  //
  // ... in a sense, we are somewhat capturing the type information within
  // this lambda, which is awesome.
  auto addConstraint =
    [descWrap, this, &instanceManifold, restricted, restrictions](size_t i){
    return [descWrap, this, &instanceManifold, restricted, restrictions, i]
    (ProblemOnManifold<T>& problem,
     const mnf::Manifold& globMani)
    {
      std::pair<typename GenericFunction<T>::intervals_t, typename Problem<T>::scaling_t>* bNSPair = &(this->boundsAndScaling_[i]);

      ::boost::shared_ptr<FunctionOnManifold<T>>
      funcOnMani = ::boost::make_shared<WrapperOnManifold<T>>
		 (*descWrap, globMani, instanceManifold, restricted, restrictions);

      while (bNSPair->first.size() < static_cast<size_t>(funcOnMani->outputSize()))
	{
	  std::cerr << "Warning: bounds not specified for output "
		    << bNSPair->first.size() << " of constrained function "
		    << funcOnMani->getName() << "." << std::endl
		    << "Assuming no bounds..." << std::endl;
	  bNSPair->first.push_back(GenericFunction<T>::makeInfiniteInterval());
	}

      while (bNSPair->second.size() < static_cast<size_t>(funcOnMani->outputSize()))
	{
	  std::cerr << "Warning: scale not specified for output "
	  << bNSPair->second.size() << " of constrained function "
	  << funcOnMani->getName() << "." << std::endl
	  << "Assuming a scale of 1..." << std::endl;
	  bNSPair->second.push_back(1.);
	}

      problem.addConstraint(funcOnMani, bNSPair->first, bNSPair->second);
    };
  }(this->boundsAndScaling_.size() - 1);

  this->lambdas_.push_back(addConstraint);

  return BoundsAndScalingSetter<T>(*bNSPair);
}

template<typename T>
template<class V, class W>
BoundsAndScalingSetter<T> ManifoldProblemFactory<T>::addConstraint(std::shared_ptr<DescriptiveWrapper<V, W>> descWrap, mnf::Manifold& instanceManifold)
{
  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;
  return this->addConstraint(descWrap,
			     instanceManifold,
			     restricted,
			     restrictions);
}

template<typename T>
BoundsAndScalingSetter<T> ManifoldProblemFactory<T>::addSum(AdderOnManifold<T>& adder)
{
  this->addElementaryManifolds(*adder.getManifold());
  this->boundsAndScaling_.push_back(std::make_pair(typename GenericFunction<T>::intervals_t(), typename Problem<T>::scaling_t()));

  std::pair<typename GenericFunction<T>::intervals_t, typename Problem<T>::scaling_t>* bNSPair = &(this->boundsAndScaling_.back());

  // We need the type of the function to instantiate the WrapperOnManifold
  // wrapper, hence why we capture the information we need here and wait
  // for the global manifold to be defined to execute this lambda.
  //
  // ... in a sense, we are somewhat capturing the type information within
  // this lambda, which is awesome.
  auto addConstraint =
    [&adder, this](size_t i){
    return [&adder, this, i]
    (ProblemOnManifold<T>& problem,
     const mnf::Manifold& globMani)
    {
      std::pair<typename GenericFunction<T>::intervals_t, typename Problem<T>::scaling_t>* bNSPair = &(this->boundsAndScaling_[i]);

      std::shared_ptr<FunctionOnManifold<T>> stdFuncOnMani = adder.getFunction(globMani);

      ROBOPTIM_DEBUG_ONLY(globMani.display();
      std::cout << "problem.function().inputSize(): " << problem.function().inputSize() << std::endl;
      std::cout << "stdFuncOnMani->inputSize(): " << stdFuncOnMani->inputSize() << std::endl;)

      ::boost::shared_ptr<FunctionOnManifold<T>>
      funcOnMani(stdFuncOnMani.get(), [stdFuncOnMani](FunctionOnManifold<T>*){});

      while (bNSPair->first.size() < static_cast<size_t>(funcOnMani->outputSize()))
	{
	  std::cerr << "Warning: bounds not specified for output "
		    << bNSPair->first.size() << " of constrained function "
		    << funcOnMani->getName() << "." << std::endl
		    << "Assuming no bounds..." << std::endl;
	  bNSPair->first.push_back(GenericFunction<T>::makeInfiniteInterval());
	}

      while (bNSPair->second.size() < static_cast<size_t>(funcOnMani->outputSize()))
	{
	  std::cerr << "Warning: scale not specified for output "
	  << bNSPair->second.size() << " of constrained function "
	  << funcOnMani->getName() << "." << std::endl
	  << "Assuming a scale of 1..." << std::endl;
	  bNSPair->second.push_back(1.);
	}

      problem.addConstraint(funcOnMani, bNSPair->first, bNSPair->second);
    };
  }(this->boundsAndScaling_.size() - 1);

  this->lambdas_.push_back(addConstraint);

  return BoundsAndScalingSetter<T>(*bNSPair);
}

// ---- //

template<typename T>
const mnf::CartesianProduct* ManifoldProblemFactory<T>::getGlobalManifold() const
{
  return constraintsManifold_.getManifold();
}

template<typename T>
std::unique_ptr<ProblemOnManifold<T>> ManifoldProblemFactory<T>::getProblem()
{
  // The global manifold is incomplete and only contains the
  // elementary manifolds of the constraints here.
  // The elementary manifolds of the objective function will
  // be added by the objLambda_ function below.
  const mnf::CartesianProduct* globalManifold = this->getGlobalManifold();

  // If no objective function was added to the problem, we add a constant
  // function to serve as as dummy one.
  if (objFunc_.numberOfFunctions() <= 0)
    {
      ROBOPTIM_DEBUG_ONLY(std::cout << "ADDING A DUMMY OBJECTIVE FUNCTION" << std::endl;)

      typename GenericConstantFunction<T>::vector_t offset (1);
      offset.setZero();

      typedef typename GenericFunction<T>::size_type size_type;
      size_type n = static_cast<size_type> (globalManifold->representationDim ());
      assert (n > 0);
      GenericConstantFunction<T>* cst = new GenericConstantFunction<T> (n, offset);

      // We make a mnf::Manifold& out of the CartesianProduct& to explicitly call
      // the overloaded constructor instead of the variadic one
      const mnf::Manifold& globberMani = *globalManifold;

      // FIXME: yet another leak to plug
      typedef DescriptiveWrapper<GenericConstantFunction<T>, ManiDesc<>> DescWrap_t;
      std::shared_ptr<DescWrap_t> descWrap = std::make_shared<DescWrap_t>(cst, globberMani);

      objFunc_.add(descWrap, globberMani);
    }

  lastObjFunc_ = objFunc_.getFunction(*globalManifold);

  std::unique_ptr<ProblemOnManifold<T>> problem
    (new ProblemOnManifold<T> (*globalManifold, toBoost(lastObjFunc_)));

  for (auto lambda : this->lambdas_)
    {
      lambda(*problem, *globalManifold);
    }

  typename GenericFunction<T>::intervals_t argBounds;

  std::function<void(const mnf::Manifold&)> setArgumentBounds =
    [this, &setArgumentBounds, &argBounds]
    (const mnf::Manifold& mani)
    {
      if (mani.isElementary())
	{
	  if (elementaryArgumentBounds_.find(mani.getInstanceId())
	      != elementaryArgumentBounds_.end())
	    {
	      auto& maniBounds = elementaryArgumentBounds_[mani.getInstanceId()];
	      argBounds.insert(argBounds.end(), maniBounds.begin(), maniBounds.end());
	    }
	  else
	    {
	      for (int i = 0; i < mani.representationDim(); ++i)
		{
		  argBounds.push_back(std::make_pair(-GenericFunction<T>::infinity(), GenericFunction<T>::infinity()));
		}
	    }
	}
      else
	{
	  for (size_t i = 0; i < mani.numberOfSubManifolds(); ++i)
	    {
	      setArgumentBounds(mani(i));
	    }
	}
    };

  setArgumentBounds(*globalManifold);

  problem->argumentBounds() = argBounds;

  return problem;
}

// ---- //

template<typename T>
template<class V, class W>
void ManifoldProblemFactory<T>::addObjective(double weight, std::shared_ptr<DescriptiveWrapper<V, W>> descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions)
{
  this->addElementaryManifolds(instanceManifold);
  objFunc_.add(weight, descWrap, instanceManifold, restricted, restrictions);
}

template<typename T>
template<class V, class W>
void ManifoldProblemFactory<T>::addObjective(double weight, std::shared_ptr<DescriptiveWrapper<V, W>> descWrap, mnf::Manifold& instanceManifold)
{
  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;

  return this->addObjective(weight,
			    descWrap,
			    instanceManifold,
			    restricted,
			    restrictions);
}

template<typename T>
template<class V, class W>
void ManifoldProblemFactory<T>::addObjective(std::shared_ptr<DescriptiveWrapper<V, W>> descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions)
{
  return this->addObjective(1.0,
			    descWrap,
			    instanceManifold,
			    restricted,
			    restrictions);
}

template<typename T>
template<class V, class W>
void ManifoldProblemFactory<T>::addObjective(std::shared_ptr<DescriptiveWrapper<V, W>> descWrap, mnf::Manifold& instanceManifold)
{
  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;

  return this->addObjective(descWrap,
			    instanceManifold,
			    restricted,
			    restrictions);
}


template<class T>
void ManifoldProblemFactory<T>::addArgumentBounds(const mnf::Manifold& manifold, const typename GenericFunction<T>::intervals_t& bounds)
{
  if (!constraintsManifold_.contains(manifold))
    {
      std::cerr << "WARNING: you are trying to add bounds on " << manifold.name() << ", a manifold not yet present in the problem's global manifold. Those bounds were still added it for later use." << std::endl;
    }

  elementaryArgumentBounds_[manifold.getInstanceId()] = bounds;

  typename GenericFunction<T>::intervals_t& storedBounds = elementaryArgumentBounds_[manifold.getInstanceId()];

  if (storedBounds.size() < static_cast<size_t>(manifold.representationDim()))
    {
      std::cerr << "WARNING: less bounds than the manifold's representation dimension were given. Completing to the correct size with (-inf, inf)." << std::endl;

      while (storedBounds.size() < static_cast<size_t>(manifold.representationDim()))
	{
	  storedBounds.push_back(std::make_pair(-GenericFunction<T>::infinity(), GenericFunction<T>::infinity()));
	}
    }

}

template<class T>
const typename GenericFunction<T>::intervals_t&
ManifoldProblemFactory<T>::getArgumentBounds(const mnf::Manifold& manifold) const
{
  // Throws if the bounds for this manifold are not found
  return elementaryArgumentBounds_.at(manifold.getInstanceId());
}

template<class T>
ManifoldProblemFactory<T>::ManifoldProblemFactory()
{
  // We only call this method here to set the objective funtion lambda,
  // thus avoiding to duplicate or isolate the code needed
  // to create it
  this->reset();
}

template<typename T>
void ManifoldProblemFactory<T>::reset()
{
  constraintsManifold_.clear();
  boundsAndScaling_.clear();
  lambdas_.clear();
  objFunc_.clear();
}

// ---- //

template<typename T>
BoundsAndScalingSetter<T>& BoundsAndScalingSetter<T>::setBounds(typename GenericFunction<T>::intervals_t& bounds)
{
  this->bNSPair_.first.clear();
  this->bNSPair_.first.insert(this->bNSPair_.first.end(), bounds.begin(), bounds.end());

  return *this;
}

template<typename T>
BoundsAndScalingSetter<T>& BoundsAndScalingSetter<T>::setScaling(typename Problem<T>::scaling_t& scaling)
{
  this->bNSPair_.second.clear();
  this->bNSPair_.second.insert(this->bNSPair_.second.end(), scaling.begin(), scaling.end());

  return *this;
}

template<typename T>
BoundsAndScalingSetter<T>::BoundsAndScalingSetter(std::pair<typename GenericFunction<T>::intervals_t, typename Problem<T>::scaling_t>& bNSPair)
  : bNSPair_(bNSPair)
{}

}

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_PROBLEM_FACTORY_HXX
