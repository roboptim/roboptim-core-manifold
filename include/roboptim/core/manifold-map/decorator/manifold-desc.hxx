#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_DESC_HXX
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_DESC_HXX

# include<manifolds/Manifold.h>
# include<manifolds/CartesianProduct.h>


template<template <typename> class T>
struct Pusher
{
  template<class FI>
  static void backPusher(std::vector<pgs::Manifold*>& manifolds, FI* function, const T<FI>& t)
  {
    // TODO: using new to create manifolds is mandatory, since we store a pointer
    // to them in the cartesian product. However, because we are not using a
    // shared pointer, those manifolds are then lost in memory leaks.
    manifolds.push_back(t.getInstance(function));
  }
};

namespace roboptim
{

  template<template <typename> class ... Typ>
  template<class U>
  pgs::Manifold* ManiDesc<Typ...>::getManifold(U* function)
  {
    std::vector<pgs::Manifold*> manifolds;

    // G++ seems to expand our arguments in the wrong order
    // when inside a call to a anonymous lambda.
    // This issue does not appear when expanding inside the
    // initialization of an array, so we use an initializer
    // list to expand the pack and pass that as an argument
    // to an inert anonymous lambda function.
    //
    // As for what it does, this is where we iterate over a
    // list of manifolds type to instantiate them, by using
    // the function passed as argument when needed.
    // The call to std::forward allows us to pass each type
    // of manifold to the pusher function separately.
    // 27 and 39 prevent errors when not passing any types.
    [] (std::initializer_list<int>)
      {}({27, (Pusher<Typ>::backPusher(manifolds, function,
	      (std::forward<Typ<U> >(Typ<U>()))), 39)...});

    if (manifolds.size() == 1)
      {
	return manifolds[0];
      }

    pgs::CartesianProduct* cartesian = new pgs::CartesianProduct();

    for (size_t i = 0; i < manifolds.size(); ++i)
      {
	cartesian->multiply(*(manifolds[i]));
      }

    return cartesian;
  }

}

#endif
