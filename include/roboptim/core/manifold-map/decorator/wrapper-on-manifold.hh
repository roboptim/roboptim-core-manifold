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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_WRAPPER_ON_MANIFOLD_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_WRAPPER_ON_MANIFOLD_HH

# include <vector>
# include <iostream>
# include <utility>
# include <type_traits>

# include <roboptim/core/detail/autopromote.hh>
# include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
# include <roboptim/core/manifold-map/decorator/function-on-manifold.hh>

# include <manifolds/Manifold.h>

# include <boost/noncopyable.hpp>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  /// \brief Maps a DescriptiveWrapper to a instance of a submanifold.
  ///
  /// \tparam T Matrix type
  template <typename T>
  class WrapperOnManifold :
    public FunctionOnManifold<T>
  {
    ROBOPTIM_DEFINE_FLAG_TYPE();
  public:
    ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_ (GenericTwiceDifferentiableFunction<T>);

    /// \brief Creates the mapping between the manifolds.
    ///
    /// Programmers should ideally not directly use this constructor, amd
    /// rely on the macros defined in the file manifold-map.hh to get the well
    /// defined WrapperOnManifold type, ready to be instantiated.
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
    explicit WrapperOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
      :
      FunctionOnManifold<T>
      (static_cast<size_type>(problemManifold.representationDim()),
       descWrap.fct().outputSize (),
       (boost::format ("%1%")
	% descWrap.fct().getName ()).str ()),
      fct_ (&descWrap.fct()),
      fctDiff_ (0),
      fctTwiceDiff_ (0),
      manifold_ (&functionManifold),
      mappingFromFunction_ ()
    {
      computeMapping(descWrap,
                     problemManifold,
                     functionManifold,
                     restrictedManifolds,
                     restrictions);
      if (fct_->template asType<GenericDifferentiableFunction<T>>())
        fctDiff_ = fct_->template castInto<GenericDifferentiableFunction<T>>();
      this->mappedGradient_ = gradient_t(static_cast<int> (this->mappingFromFunction_.size()));
      this->mappedJacobian_ = jacobian_t(descWrap.fct().outputSize(), static_cast<int> (this->mappingFromFunction_.size()));
      this->tangentMappedJacobian = Eigen::MatrixXd::Zero(descWrap.fct().outputSize(), this->tangentMappingFromFunctionSize_);
      this->mappedGradient_.setZero();
      this->mappedJacobian_.setZero();
      if (fct_->template asType<GenericTwiceDifferentiableFunction<T>>())
        fctTwiceDiff_ = fct_->template castInto<GenericTwiceDifferentiableFunction<T>>();
      this->mappedHessian_ = hessian_t(descWrap.fct().inputSize(), static_cast<int> (this->mappingFromFunction_.size()));
      this->mappedHessian_.setZero();
    }


    /// \brief Creates the mapping between the manifolds without restrictions.
    ///
    /// Programmers should ideally not directly use this constructor, amd
    /// rely on the macros defined in the file manifold-map.hh to get the well
    /// defined WrapperOnManifold type, ready to be instantiated.
    ///
    /// \tparam V input function type templating the DescriptiveWrapper.
    /// \tparam W descriptive manifold type templating the DescriptiveWrapper.
    ///
    /// \param fct instance of the DescriptiveWrapper
    /// \param problemManifold the manifold describing the whole variable vector.
    /// \param functionManifold the manifold describing the function's input vector.
    template<typename V, typename W>
    explicit WrapperOnManifold (DescriptiveWrapper<V, W>& fct,
                                 const mnf::Manifold& problemManifold,
                                 const mnf::Manifold& functionManifold)
      : WrapperOnManifold(fct, problemManifold, functionManifold,
			   std::vector<const mnf::Manifold*>(), std::vector<std::pair<long, long>>())
    {
    }

    /// \brief WrapperOnManifold destructor
    ~WrapperOnManifold();

    /// \brief Traits type.
    typedef typename parent_t::traits_t traits_t;

    /// \brief Print method.
    /// \param o output stream.
    std::ostream& print(std::ostream& o) const;

    /// \brief apply the jacobian on the manifold's tangent space
    virtual void manifold_jacobian (mnf::RefMat jacobian,
				    const_argument_ref arg)
      const;

  protected:
    void impl_compute (result_ref result, const_argument_ref x)
      const;

    template <typename V, typename W>
    void computeMapping(DescriptiveWrapper<V, W>& descWrap,
			const mnf::Manifold& problemManifold,
			const mnf::Manifold& functionManifold,
			std::vector<const mnf::Manifold*> restrictedManifolds,
			std::vector<std::pair<long, long>> restrictions);

    /// \brief map the input to the restricted problem
    ///
    /// \param argument the argument to map
    void mapArgument(const_argument_ref argument)
      const;

    /// \brief the function.
    const GenericFunction<T>* fct_;

    /// \brief the function, differentiable access.
    const GenericDifferentiableFunction<T>* fctDiff_;

    /// \brief the function, twice differentiable access.
    const GenericTwiceDifferentiableFunction<T>* fctTwiceDiff_;
    /// \brief the problem manifold.
    const mnf::Manifold* manifold_;

    /// \brief vector representing the restricted mapping
    std::vector<size_t> mappingFromFunction_;

    /// \brief array representing the restricted mapping for the tangent problem
    size_t* tangentMappingFromFunction_;
    /// \brief size of the array
    long tangentMappingFromFunctionSize_;

    /// \brief new input mapped to the restricted problem
    mutable vector_t mappedInput_;

    void impl_gradient (gradient_ref gradient,
                        const_argument_ref argument,
                        size_type functionId = 0)
      const;
    void impl_jacobian (jacobian_ref jacobian,
                        const_argument_ref arg)
      const;

    /// \brief gets the gradient from the restricted problem
    ///
    /// \param gradient the output gradient
    void unmapGradient(gradient_ref) const;

    /// \brief unmap the jacobian from the restricted problem
    ///
    /// \param jacobian the jacobian to unmap
    void unmapTangentJacobian(mnf::RefMat jacobian)
      const;

    /// \brief new gradient mapped to the restricted problem
    mutable gradient_t mappedGradient_;
    /// \brief new jacobian mapped to the restricted problem
    mutable jacobian_t mappedJacobian_;

    // Dirty copy ONLY
    mutable Eigen::MatrixXd tangentMappedJacobian;

    void impl_hessian(hessian_ref, const_argument_ref, size_type) const;

    /// \brief gets the hessian from the restricted problem
    ///
    /// \param hessian the output hessian
    void unmapHessian(hessian_ref hessian) const;

    /// \brief new hessian mapped to the restricted problem
    mutable hessian_t mappedHessian_;

    /// \brief gets the jacobian from the restricted problem
    ///
    /// \param instance the WrapperOnManifold
    /// \param jacobian the output jacobian
    void unmapJacobian(jacobian_ref jacobian) const;

    /// \brief sets the jacobian on the manifold's tangent space with the
    /// computed value
    ///
    /// \param instance the WrapperOnManifold
    inline void applyDiff() const;
  public:
    /// \brief Gets the manifold
    virtual const mnf::Manifold* getManifold() const;

    static const flag_t flags = ROBOPTIM_IS_ON_MANIFOLD;

    virtual flag_t getFlags() const;
  };

  template <typename T>
  std::ostream&
  operator<<(std::ostream& o, const WrapperOnManifold<T>& instWrap)
  {
    return instWrap.print(o);
  }

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/manifold-map/decorator/wrapper-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_WRAPPER_ON_MANIFOLD_HH
