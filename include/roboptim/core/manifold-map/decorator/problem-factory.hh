#ifndef ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HH
# define ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HH

# include <boost/shared_ptr.hpp>
# include <boost/any.hpp>

# include <roboptim/core/problem.hh>
# include <roboptim/core/manifold-map/decorator/manifold-map.hh>
# include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>

# include <manifolds/Manifold.h>

# include <map>
# include <vector>
# include <functional>

namespace roboptim
{

template<class U>
class ProblemFactory {
public:
  template<class V, class W>
  void addConstraint(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename V::intervals_t& bounds, typename U::scales_t& scales);
    template<class V, class W>
  void addConstraint(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename V::intervals_t& bounds);
    template<class V, class W>
  void addConstraint(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold, typename U::scales_t& scales);
    template<class V, class W>
  void addConstraint(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold);

  template<class V, class W>
  roboptim::ProblemOnManifold<U>* getProblem(roboptim::DescriptiveWrapper<V, W>& descWrap, pgs::Manifold& instanceManifold);

private:
  std::map<long, const pgs::Manifold*> elementaryInstanceManifolds_;
  std::vector<std::function<void(roboptim::ProblemOnManifold<U>&, const pgs::Manifold&)>> lambdas_;

  void addElementaryManifolds(const pgs::Manifold& instanceManifold);
  pgs::Manifold* getGlobalManifold();
};

#include <roboptim/core/manifold-map/decorator/problem-factory.hxx>

}

#endif //! ROBOPTIM_CORE_PLUGIN_PGSOLVER_PROBLEM_FACTORY_HH
