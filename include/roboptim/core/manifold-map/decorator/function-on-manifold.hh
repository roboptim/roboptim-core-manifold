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
# include <roboptim/core/twice-differentiable-function.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>

# include <manifolds/Manifold.h>

# include <boost/noncopyable.hpp>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  /// \brief OnManifold type checking flag
  const unsigned int ROBOPTIM_IS_ON_MANIFOLD = 1 << 16;

  /// \brief Maps a DescriptiveWrapper to a instance of a submanifold.
  ///
  /// \tparam T Matrix type
  template <typename T>
  class FunctionOnManifold :
    private boost::noncopyable,
    public GenericTwiceDifferentiableFunction<T>
  {
    ROBOPTIM_DEFINE_FLAG_TYPE();
  public:
    ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (GenericTwiceDifferentiableFunction<T>);

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
    explicit FunctionOnManifold
    (size_type inputSize, size_type outputSize, std::string name)
      :
      GenericTwiceDifferentiableFunction<T>
      (inputSize, outputSize, name)
    {}

    /// \brief FunctionOnManifold destructor
    ~FunctionOnManifold();

    /// \brief Traits type.
    typedef typename parent_t::traits_t traits_t;

    /// \brief Print method.
    /// \param o output stream.
    virtual std::ostream& print(std::ostream& o) const = 0;

    /// \brief apply the jacobian on the manifold's tangent space
    virtual void manifold_jacobian (mnf::RefMat jacobian,
				    const_argument_ref arg)
      const = 0;

  public:
    /// \brief Gets the manifold
    virtual const mnf::Manifold* getManifold() const = 0;

    static const flag_t flags = ROBOPTIM_IS_ON_MANIFOLD;

    virtual flag_t getFlags() const = 0;
  };

  template <typename T>
  std::ostream&
  operator<<(std::ostream& o, const FunctionOnManifold<T>& instWrap)
  {
    return instWrap.print(o);
  }

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/manifold-map/decorator/function-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
