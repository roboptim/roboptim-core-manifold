#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HXX

namespace roboptim
{

  template<template <typename> class ... Types>
  template<class U>
  Manifold* ManiDesc<Types...>::getManifold(U* function)
  {
    std::vector<pgs::Manifold*> manifolds;

    auto doAction = [&manifolds, &function](Types<FI>... manis)
      {
        int dummy[]{0, (Pusher<Types>::pushBack(manifolds, function, std::forward<Types<FI>>(manis)), 0)...};
      };//*/


    doAction((std::forward<Types<U> >(Types<U>()))...);

    std::stringstream sBuffer;

    for (int i = 0; i < manifolds.size(); ++i)
      {
        sBuffer << (i>0?" X ":"") + manifolds[i]->getName();
      }

    // FIXME: Here, simply return a Cartesian product of every manifold
    // We do not need a tree becase we will get rid of the structure
    // before comparing the two manifolds
    return new pgs::Manifold(sBuffer.str());
  }

}

#endif
