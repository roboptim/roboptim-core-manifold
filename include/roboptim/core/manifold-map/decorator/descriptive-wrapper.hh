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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
# include <vector>
# include <ostream>
# include <memory>
# include <boost/shared_ptr.hpp>
# include <boost/make_shared.hpp>

# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/manifold/deprecated.hh>

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
    typedef U wrappedFunction_t;

    /// \brief binds directly the function to the manifold
    ///
    /// This constructor does not use the existing descriptive solution.
    /// It should therefore only be used if it is necessary, as it can easily
    /// bring errors.
    ///
    /// \param f input function.
    /// \param m the manifold describing the function's input vector.
    DescriptiveWrapper (U* f, const mnf::Manifold& m);
    DescriptiveWrapper (const U* f, const mnf::Manifold& m);
    DescriptiveWrapper (std::shared_ptr<const U> f, std::shared_ptr<const mnf::Manifold> m);
    DescriptiveWrapper (boost::shared_ptr<const U> f, boost::shared_ptr<const mnf::Manifold> m);

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

    ~DescriptiveWrapper ();

    /// \brief const descriptive manifold getter
    const mnf::Manifold& manifold () const
    {
      return *manifold_;
    }

    /// \brief const wrapped function getter
    const U& fct ()
    {
      return *fct_;
    }

    static std::shared_ptr<DescriptiveWrapper<U, V>>
    makeUNCHECKEDDescriptiveWrapper(
        const U* fct,
        const mnf::Manifold& manifold) /*ROBOPTIM_CORE_MANIFOLD_DEPRECATED*/;

    static std::shared_ptr<DescriptiveWrapper<U, V>>
    makeUNCHECKEDDescriptiveWrapper(
        std::shared_ptr<const U> fct,
        std::shared_ptr<const mnf::Manifold> manifold);

    static std::shared_ptr<DescriptiveWrapper<U, V>>
    makeUNCHECKEDDescriptiveWrapper(
        boost::shared_ptr<const U> fct,
        boost::shared_ptr<const mnf::Manifold> manifold);

  private:
    /// \brief dimension check between the function and the descriptive manifold
    void checkDimension();

    ///\brief the function
    std::shared_ptr<const U> fct_;
    ///\brief the manifold
    std::shared_ptr<const mnf::Manifold> manifold_;
  };

  //template <typename U, typename V>
  //boost::shared_ptr<DescriptiveWrapper<U, V> >
  //descriptivewrapper (U* fct,
					//const mnf::Manifold& manifold)
  //{
    //return boost::make_shared<DescriptiveWrapper<U, V> > (fct, manifold);
  //}

  template <typename U, typename V>
  std::ostream&
  operator<<(std::ostream& o, DescriptiveWrapper<U, V>& descWrap)
  {
    o <<
      "DescriptiveWrapper's function:" << incindent <<
      iendl << "Name: " << descWrap.fct().getName() <<
      iendl << "Input size: " << descWrap.fct().inputSize() <<
      iendl << "Output size: " << descWrap.fct().outputSize() << decindent <<
      iendl << "DescriptiveWrapper's manifold:" << incindent <<
      iendl << "Name: " << descWrap.manifold().name() <<
      iendl << "Elementary: " <<
      (descWrap.manifold().isElementary() ? "yes" : "no" ) <<
      iendl << "Dimension: " << descWrap.manifold().representationDim()
				       << decindent;
    return o;
  }

  /// @}
} // end of namespace roboptim.

# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_DESCRIPTIVE_WRAPPER_HH
