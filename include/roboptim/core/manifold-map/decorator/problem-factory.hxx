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
void ProblemFactory<U>::addElementaryManifolds(const pgs::Manifold& instanceManifold)
{
  std::function<void(const pgs::Manifold&)> addElementaries =
    [this, &addElementaries]
    (const pgs::Manifold& manifold)
    {
      if (manifold.isElementary())
	{
	  // Only add the manifold if there was not any other one
	  // with the same id in the map
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
void ProblemFactory<U>::addConstraint(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename V::intervals_t& bounds, typename U::scales_t& scales)
{
  this->addElementaryManifolds(instanceManifold);

  // We need the type of the function to instantiate the FunctionOnManifold
  // wrapper, hence why we capture the information we need here and wait
  // for the global manifold to be defined to execute this lambda.
  //
  // ... in a sense, we are somewhat capturing the type information within
  // this lambda, which is awesome.
  auto addConstraint =
    [&descWrap, this, bounds, scales, &instanceManifold]
    (ProblemOnManifold<U>& problem,
     const pgs::Manifold& globMani)
    {
      ::boost::shared_ptr<FunctionOnManifold<typename V::parent_t>>
      funcOnMani(new FunctionOnManifold<typename V::parent_t>
		 (descWrap, globMani, instanceManifold)
		 );

       problem.addConstraint(funcOnMani, bounds, scales);
    };

  this->lambdas_.push_back(addConstraint);
}

template<class U>
template<class V, class W>
void ProblemFactory<U>::addConstraint(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename V::intervals_t& bounds)
{
  typename U::scales_t scales;

  for(int j = 0; j < descWrap.fct().outputSize(); ++j)
    {
      scales.push_back(1.);
    }

  this->addConstraint(descWrap, instanceManifold, bounds, scales);
}

template<class U>
template<class V, class W>
void ProblemFactory<U>::addConstraint(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename U::scales_t& scales)
{
  typename V::intervals_t bounds;

  for(int j = 0; j < descWrap.fct().outputSize(); ++j)
    {
      bounds.push_back(Function::makeInfiniteInterval());
    }

  this->addConstraint(descWrap, instanceManifold, bounds, scales);
}

template<class U>
template<class V, class W>
void ProblemFactory<U>::addConstraint(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold)
{
  typename V::intervals_t bounds;
  typename U::scales_t scales;

  for(int j = 0; j < descWrap.fct().outputSize(); ++j)
    {
      bounds.push_back(Function::makeInfiniteInterval());
      scales.push_back(1.);
    }

  this->addConstraint(descWrap, instanceManifold, bounds, scales);
}

// ---- //

template<class U>
pgs::CartesianProduct* ProblemFactory<U>::getGlobalManifold()
{
  pgs::CartesianProduct* globalManifold = new pgs::CartesianProduct();

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
  pgs::CartesianProduct* globalManifold = this->getGlobalManifold();

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
void ProblemFactory<U>::setObjective(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold)
{
  this->objLambda_ = [&descWrap, &instanceManifold, this](pgs::CartesianProduct& globMani)
    {
      std::vector<long> manifoldIds;
      std::function<void(const pgs::Manifold&)> addElementaries =
      [this, &addElementaries, &globMani, &manifoldIds]
      (const pgs::Manifold& manifold)
      {
	if (manifold.isElementary())
	  {
	    // TODO: in-lambda list to avoid multiple inclusion from the objective itself
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
      (descWrap, globMani, instanceManifold);

      return new ProblemOnManifold<U>(globMani, *objOnMani);
    };
}

template<class U>
ProblemFactory<U>::ProblemFactory()
{
  this->objLambda_ = [](pgs::CartesianProduct& globMani)
    {
      typename GenericConstantFunction<EigenMatrixDense>::vector_t offset (globMani.representationDim());
      GenericConstantFunction<EigenMatrixDense>* cst  = new GenericConstantFunction<EigenMatrixDense>(offset);

      // We make a pgs::Manifold& out of the CartesianProduct& to explicitly call
      // the overloaded constructor instead of the variadic one
      pgs::Manifold& globberMani = globMani;

      DescriptiveWrapper<GenericConstantFunction<EigenMatrixDense>, ManiDesc<>>
      descWrap(cst, globberMani);

      FunctionOnManifold<typename GenericConstantFunction<EigenMatrixDense>::parent_t>*
      objOnMani = new FunctionOnManifold<typename GenericConstantFunction<EigenMatrixDense>::parent_t>
      (descWrap, globMani, globMani);

      return new ProblemOnManifold<U>(globMani, *objOnMani);
      };
}

#endif //! ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HXX
