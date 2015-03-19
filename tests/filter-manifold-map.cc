// Copyright (C) 2013 by Thomas Moulard, AIST, CNRS, INRIA.
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

#include <roboptim/core/io.hh>
#include <roboptim/core/differentiable-function.hh>
#include <roboptim/core/util.hh>

#include <roboptim/core/function/cos.hh>
#include <roboptim/core/filter/manifold-map/instance-wrapper.hh>
#include <roboptim/core/filter/manifold-map/descriptive-wrapper.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/CartesianPower.h>
#include <manifolds/ExpMapMatrix.h>

using namespace roboptim;

struct F : public DifferentiableFunction
{
  F () : DifferentiableFunction (1, 10, "f_n (x) = n * x")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < outputSize (); ++i)
      res[i] = (value_type)i * argument[0];
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
		      size_type functionId) const
  {
    grad.setZero ();
    grad[0] = (value_type)functionId;
  }
};

boost::shared_ptr<boost::test_tools::output_test_stream> output;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (manifold_map_test)
{
  // FIXME: A test should be written here

  output = retrievePattern("filter-manifold-map");

  boost::shared_ptr<F> f (new F());

  pgs::RealSpace pos(3); pos.name() = "position";
  pgs::SO3<pgs::ExpMapMatrix> ori; ori.name() = "orientation";
  pgs::CartesianProduct freeFlyer(pos, ori);
  pgs::RealSpace joints(10); joints.name() = "joints";
  pgs::CartesianProduct robot(freeFlyer, joints);

  pgs::CartesianProduct cartProd(joints, ori);
  pgs::CartesianProduct myFuncManifold(cartProd, pos);

  InstanceWrapper<DifferentiableFunction> instWrap(f, robot, myFuncManifold);

  (*output) << instWrap;
  std::cout << instWrap;
  DescriptiveWrapper<DifferentiableFunction> descWrap(f, myFuncManifold);
  std::cout << descWrap;

  BOOST_CHECK (output->match_pattern());
}

BOOST_AUTO_TEST_SUITE_END ()
