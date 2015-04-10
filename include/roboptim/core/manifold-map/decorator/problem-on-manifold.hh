#ifndef ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_ON_MANIFOLD_HH

# include <boost/shared_ptr.hpp>

# include <roboptim/core/problem.hh>
# include <pgsolver/solver/Problem.h>
# include <roboptim/core/manifold-map/decorator/manifold-map.hh>

# include <manifolds/Manifold.h>
# include <manifolds/RealSpace.h>

namespace roboptim {

    template<class T>
    class ProblemOnManifold : public T
    {
    public:
      typedef T wrappedType_t;

      template<typename ... Types>
      ProblemOnManifold(pgs::Manifold& manifold, Types& ... args);

      pgs::Manifold& getManifold() const;

      virtual ~ProblemOnManifold();

    private:
      pgs::Manifold& manifold_;
    };

# include <roboptim/core/manifold-map/decorator/problem-on-manifold.hxx>
}

#endif //! ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_ON_MANIFOLD_HH
