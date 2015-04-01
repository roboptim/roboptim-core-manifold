#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_MAP_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_MAP_HH

# include <roboptim/core/filter/manifold-map/descriptive-wrapper.hh>
# include <roboptim/core/filter/manifold-map/manifold-desc.hh>
# include <roboptim/core/filter/manifold-map/instance-wrapper.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>

#define ROBOPTIM_DESCRIPTIVE_FWD_DECS(NAME) typedef NAME type

#define REAL_SPACE(num) Real<num>::Space
#define DESC_MANIFOLD(name, ...) typedef ManiDesc< __VA_ARGS__> name
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

#define BIND_FUNCTION_ON_MANIFOLD(function, manifold) typedef DescriptiveWrapper<function, manifold> function##_On_##manifold; \
  typedef InstanceWrapper<typename function::parent_t> Instance_##function##_On_##manifold
#define NAMED_FUNCTION_BINDING(name, function, manifold) typedef DescriptiveWrapper<function, manifold> name; \
  typedef InstanceWrapper<typename function::parent_t> Instance_##name

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

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_MAP_HH
