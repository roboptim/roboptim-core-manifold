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

  /// \brief Factory to help with the generation of the global manifold
  ///
  /// \tparam U robptim problem type to be generated.
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
    void addConstraint(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename V::intervals_t& bounds, typename U::scales_t& scales);

    template<class V, class W>
    void addConstraint(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename V::intervals_t& bounds);
    template<class V, class W>
    void addConstraint(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename U::scales_t& scales);
    template<class V, class W>
    void addConstraint(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold);

    /// \brief set the objective function of the problem
    ///
    /// \param descWrap DescriptiveWrapper of the objective function
    /// \param instanceManifold manifold on which the objective function will be evaluated
    template<class V, class W>
    void setObjective(DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold);

    /// \brief generate and return the described problem
    ProblemOnManifold<U>* getProblem();

  private:
    /// \brief elementary manifolds composing the global manifold of the problem
    std::map<long, const pgs::Manifold*> elementaryInstanceManifolds_;
    /// \brief each std::function instantiate and add a constraint to the problem
    std::vector<std::function<void(ProblemOnManifold<U>&, const pgs::Manifold&)>> lambdas_;
    /// \brief instantiate and add the objective function to the problem
    std::function<ProblemOnManifold<U>*(pgs::CartesianProduct&)> objLambda_;

    /// \brief break down a manifold in elementary manifolds and adds them to the map
    void addElementaryManifolds(const pgs::Manifold& instanceManifold);
    /// \brief assemble all elemntary manifolds into the global manifold of the problem
    pgs::CartesianProduct* getGlobalManifold();
  };

#include <roboptim/core/manifold-map/decorator/problem-factory.hxx>

  /// @}

} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HH
