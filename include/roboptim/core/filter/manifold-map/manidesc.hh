#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HH

#define ROBOPTIM_DESCRIPTIVE_FORWARD_DECS(NAME) typedef typename NAME##::parent_t type

#define REAL_SPACE(num) Real<num>::Space
#define DESC_MANIFOLD(name, ...) typedef ManiDesc< __VA_ARGS__> name
#define DEFINE_MANIFOLD(name) template<class U>\
  struct Manifold_##name{\
  static Manifold* getInstance(U* function);\
  };\
  template<class U>\
  Manifold* Manifold_##name <U>::getInstance(U* function)

namespace roboptim
{
  template<template <typename> class ... Types>
  class ManiDesc
  {
  public:

    template<class U>
    static Manifold* getManifold(U* function = nullptr);

    ROBOPTIM_DESCRIPTIVE_FORWARD_DECS(U);

  };

}
#endif
