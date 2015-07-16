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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_SUM_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_SUM_ON_MANIFOLD_HH

# include <vector>
# include <iostream>
# include <utility>
# include <type_traits>
# include <memory>
# include <functional>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
# include <roboptim/core/manifold-map/decorator/function-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/manifold-merger.hh>

# include <manifolds/Manifold.h>

# include <boost/noncopyable.hpp>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  template<typename T>
  class AdderOnManifold;

  template<typename T>
  class SumOnManifold
    : public GenericTwiceDifferentiableFunction<T>
  {
    ROBOPTIM_DEFINE_FLAG_TYPE();
  public:
    ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (GenericTwiceDifferentiableFunction<T>);

    friend class AdderOnManifold<T>;

  private:
    std::vector<std::shared_ptr<FunctionOnManifold<T> > > functions_;
    std::vector<double> weights_;

    flag_t flags_;

    mutable result_t resultBuffer;
    mutable jacobian_t jacobianBuffer;
    mutable hessian_t hessianBuffer;

    SumOnManifold(std::vector<std::shared_ptr<FunctionOnManifold<T> > >& functions,
		  std::vector<double> weights,
		  size_type inputSize, size_type outputSize, std::string name)
      : GenericTwiceDifferentiableFunction<T>(inputSize, outputSize, name),
	functions_(functions),
	weights_(weights),
	resultBuffer(outputSize),
	jacobianBuffer(inputSize, outputSize)
    {
      flags_ = functions_[0]->getFlags();

      for (auto func : functions_)
	{
	  flags_ &= func->getFlags();
	}
    }

  public:
    void impl_compute (result_ref result, const_argument_ref x)
      const
    {
      size_t i = 0;
      for (auto function : functions_)
	{
	  resultBuffer.setZero();
	  (*function)(resultBuffer, x);
	  result += weights_[i++] * resultBuffer;
	}
    }

    void impl_gradient (gradient_ref,
                        const_argument_ref,
                        size_type)
      const
    {
      std::cerr << "UNIMPLEMENTED GRADIENT IN SUM ON MANIFOLD" << std::endl;
    }

    void impl_jacobian (jacobian_ref jacobian,
                        const_argument_ref arg)
      const
    {
      size_t i = 0;
      for (auto function : functions_)
	{
	  jacobianBuffer.setZero();
	  function->jacobian(jacobianBuffer, arg);
	  jacobian += weights_[i++] * jacobianBuffer;
	}
    }

    void impl_hessian(hessian_ref hessian,
		      const_argument_ref arg,
		      size_type functionId = 0)
      const
    {
      size_t i = 0;
      for (auto function : functions_)
	{
	  hessianBuffer.setZero();
	  function->hessian(hessianBuffer, arg, functionId);
	  hessian += weights_[i++] * hessianBuffer;
	}
    }

    virtual flag_t getFlags() const
    {
      return flags_;
    }
  };

  template<typename T>
  class AdderOnManifold
  {
    typedef typename std::function<std::shared_ptr<FunctionOnManifold<T>>(const mnf::Manifold&)> descWrap_storage_t;

    std::vector<descWrap_storage_t> functionsToSum_;
    std::vector<double> weights_;
    ManifoldMerger merger_;

  public:
    AdderOnManifold(){}

    template<typename U, typename V>
    void add(DescriptiveWrapper<U, V>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions)
    {
      add(1.0, descWrap, instanceManifold, restricted, restrictions);
    }

    template<typename U, typename V>
    void add(DescriptiveWrapper<U, V>& descWrap, mnf::Manifold& instanceManifold)
    {
      std::vector<const mnf::Manifold*> restricted;
      std::vector<std::pair<long, long>> restrictions;

      add(descWrap, instanceManifold, restricted, restrictions);
    }

    template<typename U, typename V>
    void add(double weight, DescriptiveWrapper<U, V>& descWrap, mnf::Manifold& instanceManifold)
    {
      std::vector<const mnf::Manifold*> restricted;
      std::vector<std::pair<long, long>> restrictions;

      add(weight, descWrap, instanceManifold, restricted, restrictions);
    }

    template<typename U, typename V>
    void add(double weight, DescriptiveWrapper<U, V>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions)
    {
      descWrap_storage_t lambda =
	[&descWrap, this, &instanceManifold, restricted, restrictions]
	(const mnf::Manifold& globMani)
	{
	  return std::shared_ptr<FunctionOnManifold<T>>
	  (new WrapperOnManifold<T>
	   (descWrap, globMani, instanceManifold, restricted, restrictions)
	   );
	};

      functionsToSum_.push_back(lambda);
      weights_.push_back(weight);
      merger_.addManifold(instanceManifold);
    }

    std::shared_ptr<FunctionOnManifold<T>> getFunction(const mnf::Manifold& globMani)
    {
      mnf::Manifold* sumManifold = merger_.getManifold();

      std::vector<std::shared_ptr<FunctionOnManifold<T> > > functions;
      std::string name;

      for (descWrap_storage_t& lambda : functionsToSum_)
	{
	  functions.push_back(lambda(*sumManifold));
	  name += (name.size() > 0 ? " ":"") + functions.back()->getName();
	}

      typedef DescriptiveWrapper<SumOnManifold<T>, ManiDesc<>> descWrap_t;

      SumOnManifold<T>* sumFunction = new SumOnManifold<T>(functions,
							   weights_,
							   functions[0]->inputSize(),
							   functions[0]->outputSize(),
							   name);

      descWrap_t* descWrap = descWrap_t::makeUNCHECKEDDescriptiveWrapper(sumFunction, *sumManifold);

      return std::shared_ptr<FunctionOnManifold<T>>(new WrapperOnManifold<T>
						    (*descWrap, globMani, *sumManifold)
						    );

      delete descWrap;
    }

    void clear()
    {
      functionsToSum_.clear();
      weights_.clear();
      merger_.clear();
    }

    mnf::Manifold* getManifold()
    {
      return merger_.getManifold();
    }

    size_t numberOfFunctions()
    {
      return functionsToSum_.size();
    }
  };

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/manifold-map/decorator/sum-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_SUM_ON_MANIFOLD_HH
