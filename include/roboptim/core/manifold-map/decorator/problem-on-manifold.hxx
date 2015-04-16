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

template<class T>
template<typename ... Types>
ProblemOnManifold<T>::ProblemOnManifold(pgs::Manifold& manifold, Types& ... args)
  : T(args...),
    manifold_(manifold)
{
}

template<class T>
pgs::Manifold& ProblemOnManifold<T>::getManifold() const
{
  return this->manifold_;
}

template<class T>
ProblemOnManifold<T>::~ProblemOnManifold()
{
}
