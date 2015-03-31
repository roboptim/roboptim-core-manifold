#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HH



namespace roboptim
{
  template<template <typename> class ... Types>
  class ManiDesc
  {
  public:
    template<class U>
    static pgs::Manifold* getManifold(U* function = nullptr);

  };

}

# include "manifold-desc.hxx"

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_DESC_HH
