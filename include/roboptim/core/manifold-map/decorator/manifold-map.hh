// Copyright (C) 2015 by Grégoire Duchemin, AIST, CNRS, EPITA
//                       Félix Darricau, AIST, CNRS, EPITA
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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MAP_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MAP_HH

# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
# include <roboptim/core/manifold-map/decorator/manifold-desc.hh>
# include <roboptim/core/manifold-map/decorator/function-on-manifold.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>

#define ROBOPTIM_DESCRIPTIVE_FWD_DECS(NAME) typedef NAME type

#define REAL_SPACE(num) roboptim::Real<num>::Space
#define DESC_MANIFOLD(name, ...) typedef roboptim::ManiDesc< __VA_ARGS__> name
#define DEFINE_MANIFOLD_FROM_FUNCTION(name) template<class U>\
  struct Manifold_##name{\
  static pgs::Manifold* getInstance(U* function);\
  };\
  template<class U>\
  pgs::Manifold* Manifold_##name <U>::getInstance(U* function)
#define DEFINE_MANIFOLD(name) template<class U>\
  struct Manifold_##name{\
  static pgs::Manifold* getInstance(U* function);\
  };\
  template<class U>\
  pgs::Manifold* Manifold_##name <U>::getInstance(U*)

#define BIND_FUNCTION_ON_MANIFOLD(function, manifold) typedef roboptim::DescriptiveWrapper<function, manifold> function##_On_##manifold; \
  typedef roboptim::FunctionOnManifold<typename function::parent_t> Instance_##function##_On_##manifold
#define NAMED_FUNCTION_BINDING(name, function, manifold) typedef roboptim::DescriptiveWrapper<function, manifold> name; \
  typedef roboptim::FunctionOnManifold<typename function::parent_t> Instance_##name

// Library-defined elementary descriptive manifolds
// I do not think we should put those in a namespace of their own,
// because the user has to use the classes' names when (s)he defines
// a manifold.
//
// FIXME: define all remaining manifolds

namespace roboptim {

  template <class FI>
  struct S2
  {
    static pgs::Manifold* getInstance(FI*)
    {
      return new pgs::S2();
    }
  };

  template <class FI>
  struct SO3
  {
    static pgs::Manifold* getInstance(FI*)
    {
      return new pgs::SO3<pgs::ExpMapMatrix>();
    }
  };

  template<int I = 1>
  struct Real
  {
    template <class FI>
    struct Space
    {
      static pgs::Manifold* getInstance(FI*)
      {
	return new pgs::RealSpace(I);
      }
    };
  };

  template<class FI>
  struct Automated_Real
  {
  public:
    static pgs::Manifold* getInstance(FI* function)
    {
      return new pgs::RealSpace(function->getSize());
    }
  };

}

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MAP_HH
