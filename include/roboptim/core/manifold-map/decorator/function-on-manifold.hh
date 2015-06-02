// Copyright (C) 2015 by Félix Darricau, AIST, CNRS, EPITA
//                       Grégoire Duchemin, AIST, CNRS, EPITA
//
// This file is part of the roboptim.
//
// roboptim is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
# include <vector>
# include <iostream>
# include <utility>
# include <type_traits>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
# include <roboptim/core/manifold-map/decorator/helper/generic-differentiable-function-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/helper/base-function-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/helper/differentiable-function-on-manifold.hh>
# include <roboptim/core/manifold-map/decorator/helper/twice-differentiable-function-on-manifold.hh>

# include <manifolds/Manifold.h>

# include <boost/noncopyable.hpp>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  // ---- //

  // TODO: Move to a separate file
  /// \brief Maps a DescriptiveWrapper to a instance of a submanifold.
  ///
  /// WARNING: A FunctionOnManifold must be templated on the underlying roboptim
  /// function type
  ///
  /// \tparam U input roboptim function type.
  template <typename U>
  class FunctionOnManifold :
    private boost::noncopyable,
    public std::conditional<std::is_base_of<roboptim::GenericTwiceDifferentiableFunction<typename U::traits_t>,
					    U>::value,
			    TwiceDifferentiableFunctionOnManifold<U>,
			    typename std::conditional<std::is_base_of<roboptim::GenericDifferentiableFunction<typename U::traits_t>,
								      U>::value,
						      DifferentiableFunctionOnManifold<U>,
						      BaseFunctionOnManifold<U> >::type >::type
  {
  public:
    ROBOPTIM_FUNCTION_FWD_TYPEDEFS_ (U);

    /// \brief Creates the mapping between the manifolds.
    ///
    /// Programmers should ideally not directly use this constructor, amd
    /// rely on the macros defined in the file manifold-map.hh to get the well
    /// defined FunctionOnManifold type, ready to be instantiated.
    ///
    /// \tparam V input function type templating the DescriptiveWrapper.
    /// \tparam W descriptive manifold type templating the DescriptiveWrapper.
    ///
    /// \param descwrap instance of the DescriptiveWrapper
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    /// \param restrictedManifolds a list of elementary Manifolds to be restricted to a part of themselves
    /// \param restrictions the restrictions applying to the selected manifolds, represented as (startingIndex, size). If a single one is given, it will apply to all restricted manifolds.
    template <typename V, typename W>
    explicit FunctionOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
      :
      std::conditional<std::is_base_of<roboptim::GenericTwiceDifferentiableFunction<typename U::traits_t>,
      U>::value,
      TwiceDifferentiableFunctionOnManifold<U>,
      typename std::conditional<std::is_base_of<roboptim::GenericDifferentiableFunction<typename U::traits_t>,
						U>::value,
				DifferentiableFunctionOnManifold<U>,
				BaseFunctionOnManifold<U> >::type >::type(descWrap, problemManifold, functionManifold, restrictedManifolds, restrictions)
    {
    }


    /// \brief Creates the mapping between the manifolds without restrictions.
    ///
    /// Programmers should ideally not directly use this constructor, amd
    /// rely on the macros defined in the file manifold-map.hh to get the well
    /// defined FunctionOnManifold type, ready to be instantiated.
    ///
    /// \tparam V input function type templating the DescriptiveWrapper.
    /// \tparam W descriptive manifold type templating the DescriptiveWrapper.
    ///
    /// \param fct instance of the DescriptiveWrapper
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    template<typename V, typename W>
    explicit FunctionOnManifold (DescriptiveWrapper<V, W>& fct,
                                 const mnf::Manifold& problemManifold,
                                 const mnf::Manifold& functionManifold)
      : FunctionOnManifold(fct, problemManifold, functionManifold,
			   std::vector<const mnf::Manifold*>(), std::vector<std::pair<long, long>>())
    {
    }

    /// \brief Traits type.
    typedef typename parent_t::traits_t traits_t;

    std::ostream& print_(std::ostream& o);
  private:
  public:
  };

  template <typename U>
  std::ostream&
  operator<<(std::ostream& o, FunctionOnManifold<U>& instWrap)
  {
    return instWrap.print_(o);
  }

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/manifold-map/decorator/function-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
