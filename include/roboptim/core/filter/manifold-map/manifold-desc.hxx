#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HXX

# include<manifolds/Manifold.h>
# include<manifolds/CartesianProduct.h>


template<template <typename> class T>
struct Pusher
{
  template<class FI>
  static void pushBack(std::vector<pgs::Manifold*>& manifolds, FI* function, const T<FI>& t)
  {
    manifolds.push_back(t.getInstance(function));
  }
};

namespace roboptim
{

  template<template <typename> class ... Types>
  template<class U>
  pgs::Manifold* ManiDesc<Types...>::getManifold(U* function)
  {
    std::vector<pgs::Manifold*> manifolds;

    [](...){}(0, (Pusher<Types>::pushBack(manifolds, function, (std::forward<Types<U> >(Types<U>()))), 0)...);

    if (manifolds.size() == 1)
      {
	return manifolds[0];
      }

    pgs::CartesianProduct* cartesian = new pgs::CartesianProduct();

    for (size_t i = 0; i < manifolds.size(); ++i)
      {
	cartesian->multiply(*(manifolds[i]));
      }

    // FIXME: Here, simply return a Cartesian product of every manifold
    // We do not need a tree becase we will get rid of the structure
    // before comparing the two manifolds
    return cartesian;
  }

}

#endif
