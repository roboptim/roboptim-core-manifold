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
#ifndef ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HH
# define ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HH

# include <boost/shared_ptr.hpp>


# include <roboptim/core/problem.hh>
# include <roboptim/core/manifold-map/decorator/manifold-map.hh>
# include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>
# include <roboptim/core/function/constant.hh>

# include <manifolds/Manifold.h>

# include <map>
# include <vector>
# include <functional>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  template<class U>
  class ProblemFactory;

  /// \brief allows to set the bounds and scaling of a constraint
  /// after defining the constraint
  template<class U>
  struct BoundsAndScalingSetter
  {
    BoundsAndScalingSetter<U>& setBounds(typename U::function_t::intervals_t& bounds);
    BoundsAndScalingSetter<U>& setScaling(typename U::scaling_t& scaling);

  private:
    std::pair<typename U::function_t::intervals_t, typename U::scaling_t>& bNSPair_;
    BoundsAndScalingSetter(std::pair<typename U::function_t::intervals_t, typename U::scaling_t>& bNSPair);
    friend ProblemFactory<U>;
  };

  /// \brief Factory to help with the generation of the global manifold
  ///
  /// \tparam U roboptim problem type to be generated.
  template<class U>
  class ProblemFactory {
  public:

    ProblemFactory();

    /// \brief add a constraint to the problem to be generated
    ///
    /// \param descWrap DescriptiveWrapper of the constraint function
    /// \param instanceManifold manifold on which the constraint will be evaluated
    /// \param bounds (optional) the constraint's bounds
    /// \param bounds (optional) scaling parameter for the constraint
    template<class V, class W>
    BoundsAndScalingSetter<U> addConstraint(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions);
    template<class V, class W>
    BoundsAndScalingSetter<U> addConstraint(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold);

    /// \brief set the objective function of the problem
    ///
    /// \param descWrap DescriptiveWrapper of the objective function
    /// \param instanceManifold manifold on which the objective function will be evaluated
    template<class V, class W>
    void setObjective(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold, std::vector<const mnf::Manifold*>& restricted, std::vector<std::pair<long, long>>& restrictions);
    template<class V, class W>
    void setObjective(DescriptiveWrapper<V, W>& descWrap, mnf::Manifold& instanceManifold);

    /// \brief generate and return the described problem
    ProblemOnManifold<U>* getProblem();

    /// \brief reset the factory to be able to create a completely new problem
    void reset();

  private:
    /// \brief elementary manifolds composing the global manifold of the problem
    std::map<long, const mnf::Manifold*> elementaryInstanceManifolds_;
    /// \brief bounds and scaling for each constraint
    std::vector<std::pair<typename U::function_t::intervals_t, typename U::scaling_t>> boundsAndScaling_;
    /// \brief each std::function instantiate and add a constraint to the problem
    std::vector<std::function<void(ProblemOnManifold<U>&, const mnf::Manifold&)>> lambdas_;
    /// \brief instantiate and add the objective function to the problem
    std::function<ProblemOnManifold<U>*(mnf::CartesianProduct&)> objLambda_;

    /// \brief break down a manifold in elementary manifolds and adds them to the map
    void addElementaryManifolds(const mnf::Manifold& instanceManifold);
    /// \brief assemble all elementary manifolds into the global manifold of the problem
    mnf::CartesianProduct* getGlobalManifold();
  };

#include <roboptim/core/manifold-map/decorator/problem-factory.hxx>

  /// @}

} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HH
