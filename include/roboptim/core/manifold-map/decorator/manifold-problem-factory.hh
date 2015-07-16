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
#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_PROBLEM_FACTORY_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_PROBLEM_FACTORY_HH

# include <boost/shared_ptr.hpp>


# include <roboptim/core/problem.hh>
# include <roboptim/core/manifold-map/decorator/manifold-map.hh>
# include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/manifold-merger.hh>
# include <roboptim/core/manifold-map/decorator/function-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/wrapper-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/sum-on-manifold.hh>
# include <roboptim/core/function/constant.hh>

# include <manifolds/Manifold.h>

# include <map>
# include <vector>
# include <functional>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  template<typename T>
  class ManifoldProblemFactory;

  /// \brief allows to set the bounds and scaling of a constraint
  /// after defining the constraint
  template<typename T>
  struct BoundsAndScalingSetter
  {
    BoundsAndScalingSetter<T>& setBounds(typename GenericFunction<T>::intervals_t& bounds);
    BoundsAndScalingSetter<T>& setScaling(typename Problem<T>::scaling_t& scaling);

  private:
    std::pair<typename GenericFunction<T>::intervals_t, typename Problem<T>::scaling_t>& bNSPair_;
    BoundsAndScalingSetter(std::pair<typename GenericFunction<T>::intervals_t, typename Problem<T>::scaling_t>& bNSPair);
    friend ManifoldProblemFactory<T>;
  };

  /// \brief Factory to help with the generation of the global manifold
  ///
  /// \tparam T Matrix type
  template<typename T>
  class ManifoldProblemFactory {
  public:

    ManifoldProblemFactory();

    /// \brief add a constraint to the problem to be generated
    ///
    /// \param descWrap DescriptiveWrapper of the constraint function
    /// \param instanceManifold manifold on which the constraint will be evaluated
    /// \param bounds (optional) the constraint's bounds
    /// \param scaling (optional) scaling parameter for the constraint
    template<class V, class W>
    BoundsAndScalingSetter<T> addConstraint(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions);
    template<class V, class W>
    BoundsAndScalingSetter<T> addConstraint(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold);

    /// \brief This function adds a SumOnManifold function
    /// from the AdderOnManifold factory used to build it.
    BoundsAndScalingSetter<T> addSum(AdderOnManifold<T>& adder);

    /// \brief set the objective function of the problem
    ///
    /// \param descWrap DescriptiveWrapper of the objective function
    /// \param instanceManifold manifold on which the objective function will be evaluated
    template<class V, class W>
    void addObjective(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions);
    template<class V, class W>
    void addObjective(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold);
    template<class V, class W>
    void addObjective(double weight, DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions);
    template<class V, class W>
    void addObjective(double weight, DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold);


    /// \brief for a given elementary manifold, add bounds to its arguments
    ///
    /// \param manifold the elementary manifold to bound
    /// \param bounds the list of bounds for each value of a point on this manifold
    void addArgumentBounds(const mnf::Manifold& manifold, const typename GenericFunction<T>::intervals_t& bounds);

    /// \brief generate and return the described problem
    ProblemOnManifold<T>* getProblem();

    /// \brief reset the factory to be able to create a completely new problem
    void reset();

  private:
    /// \brief contains the manifolds from the constraints
    ManifoldMerger constraintsManifold_;
    /// \brief elementary manifolds composing the global manifold of the problem
    //std::set<const mnf::Manifold*> elementaryInstanceManifolds_;
    /// \brief argument bounds for each elementary manifold
    std::map<long, typename GenericFunction<T>::intervals_t> elementaryArgumentBounds_;
    /// \brief bounds and scaling for each constraint
    std::vector<std::pair<typename GenericFunction<T>::intervals_t, typename Problem<T>::scaling_t>> boundsAndScaling_;
    /// \brief each std::function instantiate and add a constraint to the problem
    std::vector<std::function<void(ProblemOnManifold<T>&, const mnf::Manifold&)>> lambdas_;
    /// \brief construct an objective function by adding functions together
    AdderOnManifold<T> objFunc_;

    /// \brief break down a manifold in elementary manifolds and adds them to the map
    void addElementaryManifolds(const mnf::Manifold& instanceManifold);
    /// \brief assemble all elementary manifolds into the global manifold of the problem
    mnf::CartesianProduct* getGlobalManifold();
  };

  /// @}

} // end of namespace roboptim

#include <roboptim/core/manifold-map/decorator/manifold-problem-factory.hxx>

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_PROBLEM_FACTORY_HH
