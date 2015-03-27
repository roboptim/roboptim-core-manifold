#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_WRAP_IN_MANIFOLD_HH
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_WRAP_IN_MANIFOLD_HH

#define ROBOPTIM_DESCRIPTIVE_FORWARD_DECS(NAME) typedef typename NAME##::parent_t type

namespace roboptim
{
  enum ManifoldType
  {
    USER  0,
    R3    1,
    S03   2
    // we should write an enum type for every useful manifold type
  }

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
