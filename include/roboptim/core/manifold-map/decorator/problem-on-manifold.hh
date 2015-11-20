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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_PROBLEM_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_PROBLEM_ON_MANIFOLD_HH

# include <memory>

# include <roboptim/core/problem.hh>
# include <roboptim/core/manifold-map/decorator/manifold-map.hh>
# include <roboptim/core/manifold/deprecated.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>

namespace roboptim
{
  /// \brief RobOptim problem on a manifold.
  ///
  /// \tparam T matrix traits.
  template<typename T>
  class ProblemOnManifold : public Problem<T>
  {
  public:
    /// \brief Constructor of a problem on manifold.
    ///
    /// \tparam Types variadic types for the problem constructor.
    /// \param manifold global manifold.
    /// \param args arguments for the constructor of Problem.
    template<typename ... Types>
    ProblemOnManifold(std::shared_ptr<const mnf::Manifold> manifold, Types&& ... args);

    /// \brief Deprecated unsafe constructor of a problem on manifold.
    template<typename ... Types>
    ProblemOnManifold(const mnf::Manifold& manifold, Types&& ... args)
    ROBOPTIM_CORE_MANIFOLD_DEPRECATED;

    virtual ~ProblemOnManifold();

    const mnf::Manifold& getManifold() const;

  protected:
    /// \brief Global manifold.
    std::shared_ptr<const mnf::Manifold> manifold_;
  };
}

# include <roboptim/core/manifold-map/decorator/problem-on-manifold.hxx>

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_PROBLEM_ON_MANIFOLD_HH
