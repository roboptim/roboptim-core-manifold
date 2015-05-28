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
#ifndef ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HXX
# define ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HXX

# include <boost/pointer_cast.hpp>



template<class U>
void ProblemFactory<U>::addElementaryManifolds(const mnf::Manifold& instanceManifold)
{
  std::function<void(const mnf::Manifold&)> addElementaries =
    [this, &addElementaries]
    (const mnf::Manifold& manifold)
    {
      if (manifold.isElementary())
	{
	  // Only add the manifold if there was not any other one
	  // with the same instanceId in the map
	  if (!this->elementaryInstanceManifolds_.count(manifold.getInstanceId()))
	    {
	      this->elementaryInstanceManifolds_[manifold.getInstanceId()] = &manifold;
	    }
	}
      else
	{
	  for (size_t i = 0; i < manifold.numberOfSubmanifolds(); ++i)
	    {
	      addElementaries(manifold(i));
	    }
	}
    };

  addElementaries(instanceManifold);
}

// ---- //

template<class U>
template<class V, class W>
BoundsAndScalesSetter<U> ProblemFactory<U>::addConstraint(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions)
{
  this->addElementaryManifolds(instanceManifold);
  this->boundsAndScales_.push_back(std::make_pair(typename V::intervals_t(), typename U::scales_t()));

  std::pair<typename U::function_t::intervals_t, typename U::scales_t>* bNSPair = &(this->boundsAndScales_.back());

  // We need the type of the function to instantiate the FunctionOnManifold
  // wrapper, hence why we capture the information we need here and wait
  // for the global manifold to be defined to execute this lambda.
  //
  // ... in a sense, we are somewhat capturing the type information within
  // this lambda, which is awesome.
  auto addConstraint =
    [&descWrap, this, &instanceManifold, restricted, restrictions](size_t i){
    return [&descWrap, this, &instanceManifold, restricted, restrictions, i]
    (ProblemOnManifold<U>& problem,
     const mnf::Manifold& globMani)
    {
      std::pair<typename U::function_t::intervals_t, typename U::scales_t>* bNSPair = &(this->boundsAndScales_[i]);

      ::boost::shared_ptr<FunctionOnManifold<typename V::parent_t>>
      funcOnMani(new FunctionOnManifold<typename V::parent_t>
		 (descWrap, globMani, instanceManifold, restricted, restrictions)
		 );

      while (bNSPair->first.size() < static_cast<size_t>(funcOnMani->outputSize()))
	{
	  std::cerr << "Warning: bounds not specified for output "
		    << bNSPair->first.size() << " of constrained function "
		    << funcOnMani->getName() << "." << std::endl
		    << "Assuming no bounds..." << std::endl;
	  bNSPair->first.push_back(U::function_t::makeInfiniteInterval());
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
  }(this->boundsAndScales_.size() - 1);

  this->lambdas_.push_back(addConstraint);

  return BoundsAndScalesSetter<U>(*bNSPair);
}

template<class U>
template<class V, class W>
BoundsAndScalesSetter<U> ProblemFactory<U>::addConstraint(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold)
{
  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;
  return this->addConstraint(descWrap,
			     instanceManifold,
			     restricted,
			     restrictions);
}

// ---- //

template<class U>
mnf::CartesianProduct* ProblemFactory<U>::getGlobalManifold()
{
  mnf::CartesianProduct* globalManifold = new mnf::CartesianProduct();

  for (auto ite = this->elementaryInstanceManifolds_.begin();
       ite != elementaryInstanceManifolds_.end();
       ++ite)
    {
      globalManifold->multiply(*(ite->second));
    }

  return globalManifold;
}

template<class U>
ProblemOnManifold<U>* ProblemFactory<U>::getProblem()
{
  // The global manifold is incomplete and only contains the
  // elementary manifolds of the cosntraints here.
  // The elementary manifolds of the objective function will
  // be added by the objLambda_ function below.
  mnf::CartesianProduct* globalManifold = this->getGlobalManifold();

  ProblemOnManifold<U>* problem = this->objLambda_(*globalManifold);

  for (auto lambda : this->lambdas_)
    {
      lambda(*problem, *globalManifold);
    }

  return problem;
}

// ---- //

template<class U>
template<class V, class W>
void ProblemFactory<U>::setObjective(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions)
{
  this->objLambda_ = [&descWrap, &instanceManifold, this, restricted, restrictions](mnf::CartesianProduct& globMani)
    {
      std::vector<long> manifoldIds;
      std::function<void(const mnf::Manifold&)> addElementaries =
      [this, &addElementaries, &globMani, &manifoldIds]
      (const mnf::Manifold& manifold)
      {
	if (manifold.isElementary())
	  {
	    if (!this->elementaryInstanceManifolds_.count(manifold.getInstanceId())
		&& std::find(manifoldIds.begin(), manifoldIds.end(), manifold.getInstanceId()) == manifoldIds.end())
	      {
		manifoldIds.push_back(manifold.getInstanceId());
		globMani.multiply(manifold);
	      }
	  }
	else
	  {
	    for (size_t i = 0; i < manifold.numberOfSubmanifolds(); ++i)
	      {
		addElementaries(manifold(i));
	      }
	  }
      };

      addElementaries(instanceManifold);

      FunctionOnManifold<typename V::parent_t>* objOnMani =
      new FunctionOnManifold<typename V::parent_t>
      (descWrap, globMani, instanceManifold, restricted, restrictions);

      return new ProblemOnManifold<U>(globMani, *objOnMani);
    };
}

template<class U>
template<class V, class W>
void ProblemFactory<U>::setObjective(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold)
{
  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;

  return this->setObjective(descWrap,
			    instanceManifold,
			    restricted,
			    restrictions);
}

template<class U>
ProblemFactory<U>::ProblemFactory()
{
  // We only call this method here to set the objective funtion lambda,
  // thus avoiding to duplicate or isolate the the code needed
  // to create it
  this->reset();
}

template<class U>
void ProblemFactory<U>::reset()
{
  elementaryInstanceManifolds_.clear();
  boundsAndScales_.clear();
  lambdas_.clear();

  this->objLambda_ = [](mnf::CartesianProduct& globMani)
    {
      typename GenericConstantFunction<typename U::function_t::traits_t>::vector_t offset (1);
      offset.setZero();
      GenericConstantFunction<typename U::function_t::traits_t>* cst  = new GenericConstantFunction<typename U::function_t::traits_t>(static_cast<typename U::function_t::size_type>(globMani.representationDim()),offset);

      // We make a mnf::Manifold& out of the CartesianProduct& to explicitly call
      // the overloaded constructor instead of the variadic one
      mnf::Manifold& globberMani = globMani;

      DescriptiveWrapper<GenericConstantFunction<typename U::function_t::traits_t>, ManiDesc<>>
      descWrap(cst, globberMani);

      FunctionOnManifold<typename GenericConstantFunction<typename U::function_t::traits_t>::parent_t>*
      objOnMani = new FunctionOnManifold<typename GenericConstantFunction<typename U::function_t::traits_t>::parent_t>
      (descWrap, globMani, globMani);

      return new ProblemOnManifold<U>(globMani, *objOnMani);
    };
}

// ---- //

template<class U>
BoundsAndScalesSetter<U>& BoundsAndScalesSetter<U>::setBounds(typename U::function_t::intervals_t& bounds)
{
  this->bNSPair_.first.clear();
  this->bNSPair_.first.insert(this->bNSPair_.first.end(), bounds.begin(), bounds.end());

  return *this;
}

template<class U>
BoundsAndScalesSetter<U>& BoundsAndScalesSetter<U>::setScales(typename U::scales_t& scales)
{
  this->bNSPair_.second.clear();
  this->bNSPair_.second.insert(this->bNSPair_.second.end(), scales.begin(), scales.end());

  return *this;
}

template<class U>
BoundsAndScalesSetter<U>::BoundsAndScalesSetter(std::pair<typename U::function_t::intervals_t, typename U::scales_t>& bNSPair)
  : bNSPair_(bNSPair)
{}

#endif //! ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HXX
