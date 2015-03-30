#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HH

#define ROBOPTIM_DESCRIPTIVE_FORWARD_DECS(NAME) typedef typename NAME##::parent_t type

#define REAL_SPACE(num) Real<num>::Space
#define DESC_MANIFOLD(name, ...) typedef ManiDesc< __VA_ARGS__> name
#define DEFINE_MANIFOLD(name) template<class FI>\
  struct Manifold_##name{\
  static Manifold* getInstance(FI* function);\
  };\
  template<class FI>\
  Manifold* Manifold_##name <FI>::getInstance(FI* function)

namespace roboptim
{
  ROBOPTIM_DESCRIPTIVE_FORWARD_DECS(U);

  template<template <typename> class ... Types>
  class ManiDesc
  {
  public:

    template<class FI>
    static Manifold* getManifold(FI* function = nullptr);

  };

  template<typename U>
  class WrapInManifold
  {
    virtual Manifold* createManifold()
    {
      // Method the programmer needs to override with its own
      return new Manifold();
    }

    template <typename ... Types>
    static U& instantiate(Types ... Args);
      // create the function of type U, using arguments Args
      // and call the createManifold method overloaded by the programmer
      // Implementation will be written in .hxx file
  }
}
