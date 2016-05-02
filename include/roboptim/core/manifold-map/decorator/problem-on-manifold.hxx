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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_PROBLEM_ON_MANIFOLD_HXX
#define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_PROBLEM_ON_MANIFOLD_HXX

# include <roboptim/core/manifold/config.hh>

namespace roboptim
{
  template<typename T>
  template<typename ... Types>
  ProblemOnManifold<T>::ProblemOnManifold
  (std::shared_ptr<const mnf::Manifold> manifold, Types&& ... args)
    : Problem<T>(args...),
      manifold_(manifold)
  {
  }

  template<typename T>
  template<typename ... Types>
  ProblemOnManifold<T>::ProblemOnManifold
  (const mnf::Manifold& manifold, Types&& ... args)
    : Problem<T>(args...),
      manifold_(&manifold, NoopDeleter<const mnf::Manifold>())
  {
  }

  template<typename T>
  ProblemOnManifold<T>::~ProblemOnManifold()
  {
  }

  template<typename T>
  const mnf::Manifold& ProblemOnManifold<T>::getManifold() const
  {
    return *manifold_;
  }

  template<typename T>
  const std::shared_ptr<const mnf::Manifold>&
  ProblemOnManifold<T>::manifold() const
  {
    return manifold_;
  }

  extern template class ROBOPTIM_CORE_MANIFOLD_DLLAPI ProblemOnManifold<EigenMatrixDense>;
  extern template class ROBOPTIM_CORE_MANIFOLD_DLLAPI ProblemOnManifold<EigenMatrixSparse>;
}
#endif //!ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_PROBLEM_ON_MANIFOLD_HXX
