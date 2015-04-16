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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MAP_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MAP_HH

# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
# include <roboptim/core/manifold-map/decorator/manifold-desc.hh>
# include <roboptim/core/manifold-map/decorator/function-on-manifold.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>
#include <manifolds/S2.h>

#define REAL_SPACE(num) roboptim::Real<num>::Space

#define DESC_MANIFOLD(name, ...) typedef roboptim::ManiDesc< __VA_ARGS__> name

#define DEFINE_MANIFOLD_FROM_FUNCTION(name) template<class U>	\
  struct Manifold_##name{					\
    static pgs::Manifold* getInstance(U* function);		\
  };								\
  template<class U>						\
  pgs::Manifold* Manifold_##name <U>::getInstance(U* function)

#define DEFINE_MANIFOLD(name) template<class U>		\
  struct Manifold_##name{				\
    static pgs::Manifold* getInstance(U* function);	\
  };							\
  template<class U>					\
  pgs::Manifold* Manifold_##name <U>::getInstance(U*)

#define BIND_FUNCTION_ON_MANIFOLD(function, manifold) typedef roboptim::DescriptiveWrapper<function, manifold> function##_On_##manifold;

#define NAMED_FUNCTION_BINDING(name, function, manifold) typedef roboptim::DescriptiveWrapper<function, manifold> name;

// Library-defined elementary descriptive manifolds
// I do not think we should put those in a namespace of their own,
// because the user has to use the classes' names when (s)he defines
// a manifold.
//
// FIXME: define all remaining manifolds


/// \brief Library-defined elementary descriptive manifolds
///
/// Each classes defined below should be instanciated by the user using one
/// of the macros defined by this file.
/// For instance, given a function type F and a manifold SO3, the user should
/// write the following sequence to define the correct descriptive manifold:
/// BIND_FUNCTION_ON_MANIFOLD(F, SO3);
///
/// He will therefore be able to define the FunctionOnManifold writing this:
/// Instance_F_On_SO3(args...) where args are the arguments needed by the
/// constructor of the FunctionOnManifold class.
///
/// Note that the classes names is the same as the manifold's one, as users
/// should never directly use them for wrapping.
namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  /// \brief Descriptive manifold for S2
  ///
  /// \tparam FI function type
  template <class FI>
  struct S2
  {
    /// \brief creates the actual manifold from the underlying library
    static pgs::Manifold* getInstance(FI*)
    {
      return new pgs::S2();
    }
  };

  /// \brief Descriptive manifold for SO3 (represented by matrices)
  ///
  /// \tparam FI function type
  template <class FI>
  struct SO3
  {
    /// \brief creates the actual manifold from the underlying library
    static pgs::Manifold* getInstance(FI*)
    {
      return new pgs::SO3<pgs::ExpMapMatrix>();
    }
  };

  /// \brief Descriptive manifold for SO3 (represented by quaternions)
  ///
  /// \tparam FI function type
  template <class FI>
  struct SO3Quat
  {
    static pgs::Manifold* getInstance(FI*)
    {
      return new pgs::SO3<pgs::ExpMapQuaternion>();
    }
  };

  /// \brief Descriptive manifold for the RealSpace family
  ///
  /// \tparam I dimension of the RealSpace
  /// \tparam FI function type
  template<int I = 1>
  struct Real
  {
    template <class FI>
    struct Space
    {
      /// \brief creates the actual manifold from the underlying library
      static pgs::Manifold* getInstance(FI*)
      {
	return new pgs::RealSpace(I);
      }
    };
  };

  /// \brief Automated descriptive manifold for the RealSpace family
  ///
  /// The dimension here is computed according to the function's size.
  ///
  /// \tparam FI function type
  template<class FI>
  struct Automated_Real
  {
    /// \brief creates the actual manifold from the underlying library
    ///
    /// \param function the function instance.
    static pgs::Manifold* getInstance(FI* function)
    {
      return new pgs::RealSpace(function->getSize());
    }
  };

  /// @}
}

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MAP_HH
