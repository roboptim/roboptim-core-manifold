// Copyright (C) 2015 by Grégoire Duchemin, AIST, CNRS, EPITA
//                       Félix Darricau, AIST, CNRS, EPITA
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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
# include <vector>
# include <ostream>
# include <boost/shared_ptr.hpp>
# include <boost/make_shared.hpp>

# include <roboptim/core/differentiable-function.hh>

# include <manifolds/Manifold.h>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  /// \brief wraps a roboptim function to a descriptive manifold
  ///
  /// \tparam U input function type.
  /// \tparam V descriptive manifold type.
  template <typename U, typename V>
  class DescriptiveWrapper
  {
  public:
    typedef boost::shared_ptr<DescriptiveWrapper> DescriptiveWrapperShPtr_t;

    /// \brief binds the function to the manifold
    ///
    /// This constructor takes a descriptive manifold as argument, that is in
    /// charge of instantiating it.
    /// It is the less error-prone way of creating a descriptive wrapper.
    ///
    /// \tparam Types list of types necessary to build the descriptive wrapper.
    ///
    /// \param args necessary arguments to build the descriptive wrapper
    template<class ... Types>
    explicit DescriptiveWrapper (Types ... args);

    /// \brief binds directly the function to the manifold
    ///
    /// This constructor does not use the existing descriptive solution.
    /// It should therefore only be used if it is necessary, as it can easily
    /// bring errors.
    ///
    /// \param f input function.
    /// \param m the manifold describing the function's input vector.
	DescriptiveWrapper (const U* f, pgs::Manifold& m);

    ~DescriptiveWrapper ();

    pgs::Manifold& manifold () const
    {
      return *manifold_;
    }

    const U& fct ()
    {
      return *fct_;
    }

  private:

    ///\brief the function
    const U*              fct_;
    ///\brief the manifold
    pgs::Manifold*  manifold_;
  };

  template <typename U, typename V>
  boost::shared_ptr<DescriptiveWrapper<U, V> >
  descriptivewrapper (U* fct,
		      pgs::Manifold& manifold)
  {
    return boost::make_shared<DescriptiveWrapper<U, V> > (fct, manifold);
  }

  template <typename U, typename V>
  std::ostream&
  operator<<(std::ostream& o, DescriptiveWrapper<U, V>& descWrap)
  {
    o <<
    "Displaying info about function in DescriptiveWrapper :" << "\n" <<
    "the function is " << descWrap.fct().getName() << "\n" <<
    "the function input size is " << descWrap.fct().inputSize() << "\n" <<
    "the function output size is " << descWrap.fct().outputSize() << "\n" <<
    "Displaying info about the function manifold in DescriptiveWrapper :" << "\n" <<
    "the manifold is " << descWrap.manifold().name() << "\n" <<
    "Is the manifold elementary ? : " <<
    (descWrap.manifold().isElementary() ? "yes" : "no" ) << "\n" <<
    "the manifold dimension is " << descWrap.manifold().representationDim() << "\n";
    return o;
  }

  /// @}
} // end of namespace roboptim.

# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
