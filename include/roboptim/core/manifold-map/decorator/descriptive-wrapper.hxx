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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HXX
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HXX
# include <boost/format.hpp>
# include <manifolds/Manifold.h>
# include "manifold-desc.hh"

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  template <typename U, typename V>
  template<class ... Types>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (Types ... args)
  {
    this->fct_ = new U(args...);
    this->manifold_ = V::getManifold(this->fct_);
    checkDimension();
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (const U* fct, pgs::Manifold& manifold)
  : fct_ (fct),
    manifold_ (&manifold)
  {
    checkDimension();
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::~DescriptiveWrapper()
  {
  }

  template <typename U, typename V>
  void DescriptiveWrapper<U, V>::checkDimension()
  {
    long size = manifold_->representationDim();
    if (fct_->inputSize() != size)
    {
      std::stringstream* error = new std::stringstream;
      (*error) << "Representation dims mismatch on manifold "
               << manifold_->name() << " using function "
               << fct_->getName() << ". Expected dimension :" << size
               << ", actual one is " << fct_->inputSize();
      throw std::runtime_error (error->str());
    }
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HXX
