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
#include <sstream>

#include <boost/format.hpp>

#include <roboptim/core/util.hh>
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
    : roboptim::GenericDifferentiableFunction<T>
      (n, 1, (boost::format ("f(x) = %1% * sum(x)") % k).str ()),
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

  using roboptim::operator<<;

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
  // TODO: use different sizes
  mnf::RealSpace f0Manif (n);
  mnf::RealSpace f1Manif (n);
  mnf::CartesianProduct pbManif;
  pbManif.multiply(f0Manif).multiply(f1Manif);

  std::stringstream ss;
  ss << "F0_R" << n;
  f0Manif.name () = ss.str ();
  ss.str ("");
  ss << "F1_R" << n;
  f1Manif.name () = ss.str ();
  ss.str ("");

  (*output) << "F0 manifold: " << f0Manif.name () << std::endl;
  (*output) << "F1 manifold: " << f1Manif.name () << std::endl;
  (*output) << "Problem manifold: " << pbManif.name () << std::endl;

  ROBOPTIM_DESC_MANIFOLD(f0Manif_t, ROBOPTIM_REAL_SPACE(n));
  ROBOPTIM_DESC_MANIFOLD(f1Manif_t, ROBOPTIM_REAL_SPACE(n));
  ROBOPTIM_NAMED_FUNCTION_BINDING(Desc_F0_Manif, testF_t, f0Manif_t);
  ROBOPTIM_NAMED_FUNCTION_BINDING(Desc_F1_Manif, testF_t, f1Manif_t);

  std::shared_ptr<Desc_F0_Manif> descWrapF0 = std::make_shared<Desc_F0_Manif>(n, k0);
  std::shared_ptr<Desc_F1_Manif> descWrapF1 = std::make_shared<Desc_F1_Manif>(n, k1);

  (*output) << *descWrapF0 << std::endl;
  (*output) << *descWrapF1 << std::endl;

  // Create Adder helper and add the functions
  typedef roboptim::AdderOnManifold<T> adder_t;
  adder_t adder;
  typename adder_t::functionPtr_t fSum;
  const value_type w0 = 0.2;
  const value_type w1 = 0.8;
  adder.add (w0, descWrapF0, f0Manif);
  adder.add (w1, descWrapF1, f1Manif);
  BOOST_CHECK (adder.numberOfFunctions () == 2);
  (*output) << "Adder manifold: "
            << adder.getManifold ()->name () << std::endl;
  (*output) << "Adder manifold dimension: "
            << adder.getManifold ()->representationDim () << std::endl;

  // Get the sum of these functions on manifolds
  fSum = adder.getFunction (pbManif);
  (*output) << *fSum << std::endl;

  // Yes, that used to fail...
  BOOST_CHECK (adder.getManifold ()->representationDim () == adder.getManifold ()->representationDim ());
  BOOST_CHECK (adder.getManifold ()->name () == adder.getManifold ()->name ());

  adder.clear();
  BOOST_CHECK (adder.numberOfFunctions () == 0);
  BOOST_CHECK_THROW (fSum = adder.getFunction (pbManif), std::runtime_error);

  // Test evaluations
  typedef roboptim::Function::matrix_t denseMatrix_t;
  typename testF_t::argument_t x0 (n);
  typename testF_t::argument_t x1 (n);
  typename testF_t::argument_t x (x0.size () + x1.size ());
  x0.setConstant (1.);
  x1.setConstant (1.);
  x << x0, x1;

  (*output) << "x0 = " << x0 << '\n';
  (*output) << "x1 = " << x1 << '\n';
  (*output) << "x = " << x << '\n';

  (*output) << "### SumOnManifold ###\n";
  (*output) << "Sum(x) = " << (*fSum) (x) << '\n';
  (*output) << "Sum.G(x) = " << denseMatrix_t (fSum->gradient (x, 0)) << '\n';
  (*output) << "Sum.J(x) = " << denseMatrix_t (fSum->jacobian (x)) << '\n';

  (*output) << "### Expected ###\n";
  (*output) << "f0(x0) = " << (*f0) (x0) << '\n';
  (*output) << "f1(x1) = " << (*f1) (x1) << '\n';
  (*output) << "(" << w0 << " f0 + " << w1 << " f1)(x) = "
            << w0 * (*f0) (x0) + w1 * (*f1) (x1) << '\n';
  (*output) << "(" << w0 << " f0).G(x0) = "
            << denseMatrix_t (w0 * f0->gradient (x0, 0)) << '\n';
  (*output) << "(" << w1 << " f1).G(x1) = "
            << denseMatrix_t (w1 * f1->gradient (x1, 0)) << '\n';
  (*output) << "(" << w0 << " f0).J(x0) = "
            << denseMatrix_t (w0 * f0->jacobian (x0)) << '\n';
  (*output) << "(" << w1 << " f1).J(x1) = "
            << denseMatrix_t (w1 * f1->jacobian (x1)) << '\n';

  std::cout << output->str() << std::endl;
  BOOST_CHECK (output->match_pattern());
}

BOOST_AUTO_TEST_SUITE_END ()
