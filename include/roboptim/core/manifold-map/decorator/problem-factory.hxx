#ifndef ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HXX
# define ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HXX

#include <boost/pointer_cast.hpp>

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
void ProblemFactory<U>::addConstraint(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename V::intervals_t& bounds, typename U::scales_t& scales)
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
    (roboptim::ProblemOnManifold<U>& problem,
     const pgs::Manifold& globMani)
    {
      ::boost::shared_ptr<roboptim::FunctionOnManifold<typename V::parent_t>>
      funcOnMani(new roboptim::FunctionOnManifold<typename V::parent_t>
		 (descWrap, globMani, instanceManifold)
		 );

       problem.addConstraint(funcOnMani, bounds, scales);
    };

  this->lambdas_.push_back(addConstraint);
}

template<class U>
template<class V, class W>
void ProblemFactory<U>::addConstraint(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename V::intervals_t& bounds)
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
void ProblemFactory<U>::addConstraint(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename U::scales_t& scales)
{
  typename V::intervals_t bounds;

  for(int j = 0; j < descWrap.fct().outputSize(); ++j)
    {
      bounds.push_back(roboptim::Function::makeInfiniteInterval());
    }

  this->addConstraint(descWrap, instanceManifold, bounds, scales);
}

template<class U>
template<class V, class W>
void ProblemFactory<U>::addConstraint(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold)
{
  typename V::intervals_t bounds;
  typename U::scales_t scales;

  for(int j = 0; j < descWrap.fct().outputSize(); ++j)
    {
      bounds.push_back(roboptim::Function::makeInfiniteInterval());
      scales.push_back(1.);
    }

  this->addConstraint(descWrap, instanceManifold, bounds, scales);
}

// ---- //

template<class U>
pgs::Manifold* ProblemFactory<U>::getGlobalManifold()
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
template<class V, class W>
roboptim::ProblemOnManifold<U>* ProblemFactory<U>::getProblem(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold)
{
  this->addElementaryManifolds(instanceManifold);

  pgs::Manifold* globalManifold = this->getGlobalManifold();

  roboptim::FunctionOnManifold<typename V::parent_t>* instObjFunc = new roboptim::FunctionOnManifold<typename V::parent_t>(descWrap, *globalManifold, instanceManifold);

  roboptim::ProblemOnManifold<U>* problem = new roboptim::ProblemOnManifold<U>(*globalManifold,
									       *instObjFunc);

  for (auto lambda : this->lambdas_)
    {
      lambda(*problem, *globalManifold);
    }

  return problem;
}

#endif //! ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HXX
