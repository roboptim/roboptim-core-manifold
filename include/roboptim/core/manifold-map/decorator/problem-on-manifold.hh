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
#ifndef ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_ON_MANIFOLD_HH

# include <boost/shared_ptr.hpp>

# include <roboptim/core/problem.hh>
# include <roboptim/core/manifold-map/decorator/manifold-map.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>

namespace roboptim {

  template<class T>
  class ProblemOnManifold : public T
  {
  public:
    typedef T wrappedType_t;

    template<typename ... Types>
    ProblemOnManifold(mnf::Manifold& manifold, Types& ... args);

    mnf::Manifold& getManifold() const;

    virtual ~ProblemOnManifold();

  private:
    mnf::Manifold& manifold_;
  };

# include <roboptim/core/manifold-map/decorator/problem-on-manifold.hxx>
}

#endif //! ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_ON_MANIFOLD_HH
