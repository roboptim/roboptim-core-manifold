// Copyright (C) 2015 by Benjamin Chrétien, CNRS-UM LIRMM.
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

#include "shared-tests/fixture.hh"

#include <iostream>
#include <memory>

#include <roboptim/core/differentiable-function.hh>

#include <roboptim/core/manifold-map/decorator/manifold-map.hh>
#include <roboptim/core/manifold-map/decorator/sum-on-manifold.hh>
#include <roboptim/core/manifold-map/decorator/wrapper-on-manifold.hh>

#include <manifolds/RealSpace.h>

typedef boost::mpl::list< ::roboptim::EigenMatrixDense,
			  ::roboptim::EigenMatrixSparse> functionTypes_t;

template<class T>
struct F : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  F (size_type n, value_type k)
    : roboptim::GenericDifferentiableFunction<T> (n, 1, "f (x) = k Σ(x_i)"),
      k_ (k)
  {}

  void impl_compute (result_ref res, const_argument_ref x) const
  {
    res[0] = k_ * x.sum ();
  }

  void impl_gradient (gradient_ref grad, const_argument_ref, size_type) const
  {
    grad.setConstant (k_);
  }

private:
  value_type k_;
};

template <>
inline void
F<roboptim::EigenMatrixSparse>::impl_gradient(gradient_ref grad,
                                              const_argument_ref,
                                              size_type) const
{
  grad.setZero ();
  for (size_type i = 0; i < inputSize (); ++i)
    {
      grad.insert (i) = k_;
    }
}

boost::shared_ptr<boost::test_tools::output_test_stream> output;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (sum_on_manifold_test_0, T, functionTypes_t)
{
  output = retrievePattern("sum-on-manifold");

  typedef F<T> testF_t;
  typedef typename testF_t::value_type value_type;

  // Create dummy functions
  static const int n = 5;
  const value_type k0 = 2.;
  const value_type k1 = 3.;
  std::shared_ptr<testF_t> f0 = std::make_shared<testF_t> (n, k0);
  std::shared_ptr<testF_t> f1 = std::make_shared<testF_t> (n, k1);

  (*output) << *f0 << std::endl;
  (*output) << *f1 << std::endl;

  // Define manifolds
  mnf::RealSpace fManif (n);
  fManif.name () = "F manifold";
  mnf::RealSpace pbManif (1);
  pbManif.name () = "Problem manifold";

  ROBOPTIM_DESC_MANIFOLD(fManif_t, ROBOPTIM_REAL_SPACE(n));
  ROBOPTIM_NAMED_FUNCTION_BINDING(Desc_F_Manif, testF_t, fManif_t);

  Desc_F_Manif descWrapF0 (n, k0);
  Desc_F_Manif descWrapF1 (n, k1);

  mnf::RealSpace dummy (n);

  (*output) << descWrapF0 << std::endl;
  (*output) << descWrapF1 << std::endl;

  // Create Adder helper and add the functions
  typedef roboptim::AdderOnManifold<T> adder_t;
  typedef F<T> testF_t;
  adder_t adder;
  typename adder_t::functionPtr_t fSum;
  BOOST_CHECK_THROW (fSum = adder.getFunction (pbManif), std::runtime_error);
  adder.add (0.2, descWrapF0, dummy);
  adder.add (0.8, descWrapF1, dummy);

  // Get the sum of these functions on manifolds
  fSum = adder.getFunction (pbManif);
  (*output) << *fSum << std::endl;

  std::cout << output->str() << std::endl;
  BOOST_CHECK (output->match_pattern());
}

BOOST_AUTO_TEST_SUITE_END ()
