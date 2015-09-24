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

  /// \brief This RobOptim function computes a weighted sum of
  /// FunctionOnManifold. It is defined on a manifold which
  /// merges all elementary manifolds of the functions it sums
  /// as merged by a ManifoldMerger.
  ///
  /// Since messing up the arguments and wrapping around it is
  /// easy, it can only be constructed through an appropriate
  /// factory, AdderOnManifold.
  ///
  /// \tparam T RobOptim Eigen matrix type
  template<typename T>
  class SumOnManifold
    : public GenericTwiceDifferentiableFunction<T>
  {
    ROBOPTIM_DEFINE_FLAG_TYPE();
  public:
    ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (GenericTwiceDifferentiableFunction<T>);

    typedef FunctionOnManifold<T> function_t;
    typedef std::shared_ptr<function_t> functionPtr_t;

    /// This function can only be created through a factory
    friend class AdderOnManifold<T>;

  private:
    /// \brief shared pointers to the functions to sum
    std::vector<functionPtr_t> functions_;
    /// \brief weight of each function in the sum
    std::vector<double> weights_;

    /// \brief the type flags of this function sum
    flag_t flags_;

    /// \brief a buffer for the result of computing this sum
    mutable result_t resultBuffer_;
    /// \brief a buffer for the result of computing the gradient of this sum
    mutable gradient_t gradientBuffer_;
    /// \brief a buffer for the result of computing the jacobian of this sum
    mutable jacobian_t jacobianBuffer_;
    /// \brief a buffer for the result of computing the hessian of this sum
    mutable hessian_t hessianBuffer_;

    /// \brief constructs a SumOnManifold
    ///
    /// \param functions a vector of functions to be summed
    /// \param weights a vector of weights, one for each function
    /// \param inputSize the inputSize of this sum (representationDim of the merged manifold)
    /// \param outputSize the outputSize of all summed functions
    /// \param name the name of this sum
    SumOnManifold(std::vector<functionPtr_t>& functions,
		  std::vector<double> weights,
		  size_type inputSize, size_type outputSize, std::string name);

  public:
    void impl_compute (result_ref result, const_argument_ref x)
      const;

    void impl_gradient (gradient_ref,
                        const_argument_ref,
                        size_type)
      const;

    void impl_jacobian (jacobian_ref jacobian,
                        const_argument_ref x)
      const;

    void impl_hessian(hessian_ref hessian,
		      const_argument_ref x,
		      size_type functionId = 0)
      const;

    virtual flag_t getFlags() const;
  };

  /// \brief This class is a factory used to create a weighted
  /// sum of FunctionOnManifold functions.
  ///
  /// \tparam T RobOptim Eigen matrix type
  template<typename T>
  class AdderOnManifold
  {
  public:
    typedef FunctionOnManifold<T> function_t;
    typedef std::shared_ptr<function_t> functionPtr_t;

  private:
    typedef typename std::function<functionPtr_t (const mnf::Manifold&)>
      descWrap_storage_t;

    /// \brief DescriptiveWrappers of functions to be summed.
    /// Stored as lambda functions to instantiate the FunctionOnManifold
    /// upon receiving the merged manifold.
    std::vector<descWrap_storage_t> functionsToSum_;
    /// \brief weights of each function
    std::vector<double> weights_;
    /// \brief utility member to compute the merged manifold of all summed functions
    ManifoldMerger merger_;

  public:
    /// \brief Default constructor
    AdderOnManifold(){}

    /// \brief Adds a function to the sum
    ///
    /// \tparam U type of the function being added
    /// \tparam V ManiDesc of the DescriptiveWrapper being added. Irrelevant here.
    ///
    /// \param weight weight fo the function in the weighted sum
    /// \param descWrap a DescriptiveWrapper<U, V> representing the function to be added to the sum
    /// \param instanceManifold the manifold on which the function will be instantiated
    /// \param restricted list of restricted manifolds
    /// \param restrictions list of restrictions for the restricted manifolds
    template<typename U, typename V>
    void add(double weight, DescriptiveWrapper<U, V>& descWrap, const mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions);
    template<typename U, typename V>
    void add(DescriptiveWrapper<U, V>& descWrap, const mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions);

    template<typename U, typename V>
    void add(DescriptiveWrapper<U, V>& descWrap, const mnf::Manifold& instanceManifold);

    template<typename U, typename V>
    void add(double weight, DescriptiveWrapper<U, V>& descWrap, const mnf::Manifold& instanceManifold);

    /// \brief Instantiate and return a shared pointer to a function computing the
    /// weighted sum of all functions added through the add(...) methods.
    ///
    /// \param globMani the manifold of the problem on which the function will be evaluated
    functionPtr_t getFunction(const mnf::Manifold& globMani) const;

    /// \brief clears out all fields of this instance
    void clear();

    /// \brief returns the merged manifold of all added functions
    const mnf::Manifold* getManifold() const;

    /// \brief return the number of functions added to the sum
    size_t numberOfFunctions();
  };

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/manifold-map/decorator/sum-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_SUM_ON_MANIFOLD_HH
