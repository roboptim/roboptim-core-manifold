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

# include <manifolds/Manifold.h>

# include <boost/noncopyable.hpp>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  /// \brief OnManifold type checking flag
  const unsigned int roboptimIsOnManifold = 1 << 16;

  /// \brief Maps a DescriptiveWrapper to a instance of a submanifold.
  ///
  /// \tparam T Matrix type
  template <typename T>
  class FunctionOnManifold :
    private boost::noncopyable,
    public GenericTwiceDifferentiableFunction<T>
  {
    ROBOPTIM_DEFINE_FLAG_TYPE()
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
    template <typename V, typename W>
    explicit FunctionOnManifold
    (DescriptiveWrapper<V, W>& descWrap,
     const mnf::Manifold& problemManifold,
     const mnf::Manifold& functionManifold,
     std::vector<const mnf::Manifold*> restrictedManifolds,
     std::vector<std::pair<long, long>> restrictions)
      :
      GenericTwiceDifferentiableFunction<T>
      (static_cast<size_type>(problemManifold.representationDim()),
       descWrap.fct().outputSize (),
       (boost::format ("%1%")
	% descWrap.fct().getName ()).str ()),
      fct_ (&descWrap.fct()),
      manifold_ (&functionManifold)
    {
      computeMapping(descWrap,
                     problemManifold,
                     functionManifold,
                     restrictedManifolds,
                     restrictions);
      this->mappedGradient_ = gradient_t(static_cast<int> (this->mappingFromFunctionSize_));
      this->mappedJacobian_ = jacobian_t(descWrap.fct().outputSize(), static_cast<int> (this->mappingFromFunctionSize_));
      this->tangentMappedJacobian = Eigen::MatrixXd::Zero(descWrap.fct().outputSize(), this->tangentMappingFromFunctionSize_);
      this->mappedGradient_.setZero();
      this->mappedJacobian_.setZero();
      this->mappedHessian_ = hessian_t(descWrap.fct().inputSize(), static_cast<int> (this->mappingFromFunctionSize_));
      this->mappedHessian_.setZero();
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

    /// \brief FunctionOnManifold destructor
    ~FunctionOnManifold();

    /// \brief Traits type.
    typedef typename parent_t::traits_t traits_t;

    std::ostream& print_(std::ostream& o);

    /// \brief apply the jacobian on the manifold's tangent space
    void manifold_jacobian (mnf::RefMat jacobian,
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
    mutable GenericFunction<T>* fct_;
    /// \brief the problem manifold.
    const mnf::Manifold* manifold_;

    /// \brief array representing the restricted mapping
    size_t* mappingFromFunction_;
    /// \brief size of the array
    long mappingFromFunctionSize_;

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
    /// \param instance the FunctionOnManifold
    /// \param jacobian the output jacobian
    void unmapJacobian(jacobian_ref jacobian) const;

    /// \brief sets the jacobian on the manifold's tangent space with the
    /// computed value
    ///
    /// \param instance the FunctionOnManifold
    void applyDiff() const;
  public:
    /// \brief Gets the manifold
    const mnf::Manifold* getManifold() const;

    static const flag_t flags = roboptimIsOnManifold;

    virtual flag_t getFlags() const;

  };

  template <typename T>
  std::ostream&
  operator<<(std::ostream& o, FunctionOnManifold<T>& instWrap)
  {
    return instWrap.print_(o);
  }

  /// @}

} // end of namespace roboptim.


# include <roboptim/core/manifold-map/decorator/function-on-manifold.hxx>
#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_FUNCTION_ON_MANIFOLD_HH
