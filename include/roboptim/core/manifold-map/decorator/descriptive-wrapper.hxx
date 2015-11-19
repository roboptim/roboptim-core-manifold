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
# include <manifolds/RealSpace.h>
# include <roboptim/core/manifold-map/decorator/manifold-desc.hh>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  template<typename U>
  struct NoopDeleter
  {
    inline void operator() (const U*) const {}
  };

  template <typename T>
  std::shared_ptr<T> make_shared_ptr(boost::shared_ptr<T>& ptr)
  {
    return std::shared_ptr<T>(ptr.get(), [ptr](T*) mutable
                              {
                                ptr.reset();
                              });
  }

  template <typename U, typename V>
  template <class... Types>
  DescriptiveWrapper<U, V>::DescriptiveWrapper(Types... args)
      : fct_(new U(args...)),
        manifold_(V::getManifold(this->fct_.get()),
                  NoopDeleter<mnf::Manifold>())
  {
    checkDimension();
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (U* fct, const mnf::Manifold& manifold)
    : fct_(fct, NoopDeleter<U>()),
      manifold_(&manifold, NoopDeleter<mnf::Manifold>())
  {
    checkDimension();
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (const U* fct, const mnf::Manifold& manifold)
    : fct_(fct, NoopDeleter<U>()),
      manifold_(&manifold, NoopDeleter<mnf::Manifold>())
  {
    checkDimension();
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (std::shared_ptr<const U> fct, std::shared_ptr<const mnf::Manifold> manifold)
    : fct_(fct),
      manifold_(manifold)
  {
    checkDimension();
  }

  template <typename U, typename V>
  DescriptiveWrapper<U, V>::DescriptiveWrapper
  (boost::shared_ptr<const U> fct, boost::shared_ptr<const mnf::Manifold> manifold)
    : fct_(make_shared_ptr<const U>(fct)),
      manifold_(make_shared_ptr<const mnf::Manifold>(manifold))
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
	std::stringstream error;
	error << "Representation dims mismatch on manifold "
	      << manifold_->name() << " using function "
	      << fct_->getName() << ". Expected dimension :" << size
	      << ", actual one is " << fct_->inputSize();
	throw std::runtime_error (error.str());
      }
  }

  template <typename U, typename V>
  std::shared_ptr<DescriptiveWrapper<U, V>> DescriptiveWrapper<U, V>::makeUNCHECKEDDescriptiveWrapper(U* fct, const mnf::Manifold& manifold)
  {
    // FIXME: remove this garbage. Anyone currently relying on this method
    // should really, really, **REALLY** take a look at this...
    mnf::RealSpace realR(fct->inputSize());
    const mnf::Manifold* m = &realR;
    std::shared_ptr<DescriptiveWrapper<U, V>> descWrap = std::make_shared<DescriptiveWrapper<U, V>>(fct, *m);

    descWrap->manifold_.reset(&manifold, NoopDeleter<mnf::Manifold>());

    return descWrap;
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HXX
