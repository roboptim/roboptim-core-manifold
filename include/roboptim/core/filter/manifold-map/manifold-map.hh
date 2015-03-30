#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_MAP_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_MAP_HH

# include <roboptim/core/filter/manifold-map/descriptive-wrapper.hh>
# include <roboptim/core/filter/manifold-map/manifold-desc.hh>
# include <roboptim/core/filter/manifold-map/instance-wrapper.hh>

#define ROBOPTIM_DESCRIPTIVE_FORWARD_DECS(NAME) typedef typename NAME##::parent_t type

#define REAL_SPACE(num) Real<num>::Space
#define DESC_MANIFOLD(name, ...) typedef ManiDesc< __VA_ARGS__> name
#define DEFINE_MANIFOLD(name) template<class U>\
  struct Manifold_##name{\
  static pgs::Manifold* getInstance(U* function);\
  };\
  template<class U>\
  pgs::Manifold* Manifold_##name <U>::getInstance(U* function)

#define BIND_FUNCTION_ON_MANIFOLD(function, manifold) typedef DW<function, manifold> function##_On_##manifold ;
#define DEFAULT_FUNCTION_BINDING(function, manifold) typedef DW<function, manifold> Wrapped##function ;

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_MANIFOLD_MAP_HH
