

template<class T>
template<typename ... Types>
ProblemOnManifold<T>::ProblemOnManifold(pgs::Manifold& manifold, Types& ... args)
  : T(args...),
    manifold_(manifold)
{
}

template<class T>
pgs::Manifold& ProblemOnManifold<T>::getManifold() const
{
  return this->manifold_;
}

template<class T>
ProblemOnManifold<T>::~ProblemOnManifold()
{
}
