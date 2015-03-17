#ifndef ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_HXX
# define ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_HXX
# include <boost/format.hpp>
# include <vector>
# include <queue>
# include <manifolds/Manifold.h>

namespace roboptim
{

  template <typename U>
  ManifoldMap<U>::ManifoldMap
  (boost::shared_ptr<U> origin,
   pgs::Manifold& problemManifold,
   pgs::Manifold& functionManifold)
    : detail::AutopromoteTrait<U>::T_type
      (origin->inputSize (),
       origin->outputSize (),
       (boost::format ("%1%")
	% origin->getName ()).str ()),
      origin_ (origin)
  {
    this->mappingSize_ = problemManifold.dim();
    this->mapping_ = new int[this->mappingSize_];

    // Retrieve the names of all manifolds to match in the problem
    std::vector<std::string> manifoldsToFind(functionManifold.numberOfSubmanifolds());

    auto retrieveManifoldNames = [&manifoldsToFind](pgs::Manifold& manifold)
      {

	manifoldsToFind.push_back(manifold.name());

	for(size_t i = 0; i < manifold.numberOfSubmanifolds(); ++i)
	  {
	    retrieveManifoldNames(manifold(i));
	  }

      };

    retrieveManifoldNames(problemManifold);

    // Traverse the problem's manifold, filling the mapping_ array


    // If manifoldsToFind is not empty, at least one submanifold could not be matched
    // in the problem.
    // TODO: raise an exception.


  }

  template <typename U>
  ManifoldMap<U>::~ManifoldMap()
  {
    delete this->mapping_;
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_FILTER_MANIFOLD_MAP_HXX
