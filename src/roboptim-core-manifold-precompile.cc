# include <roboptim/core/manifold-map/decorator/function-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/manifold-problem-factory.hh>
# include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/sum-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/wrapper-on-manifold.hh>

namespace roboptim
{

  template class FunctionOnManifold<EigenMatrixDense>;
  template class FunctionOnManifold<EigenMatrixSparse>;

  template struct BoundsAndScalingSetter<EigenMatrixDense>;
  template struct BoundsAndScalingSetter<EigenMatrixSparse>;

  template class ManifoldProblemFactory<EigenMatrixDense>;
  template class ManifoldProblemFactory<EigenMatrixSparse>;

  template class ProblemOnManifold<EigenMatrixDense>;
  template class ProblemOnManifold<EigenMatrixSparse>;

  template class SumOnManifold<EigenMatrixDense>;
  template class SumOnManifold<EigenMatrixSparse>;

  template class AdderOnManifold<EigenMatrixDense>;
  template class AdderOnManifold<EigenMatrixSparse>;

  template class WrapperOnManifold<EigenMatrixDense>;
  template class WrapperOnManifold<EigenMatrixSparse>;
}
