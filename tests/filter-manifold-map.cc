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

#include <roboptim/core/differentiable-function.hh>

#include <roboptim/core/filter/manifold-map/instance-wrapper.hh>
#include <roboptim/core/filter/manifold-map/descriptive-wrapper.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>

using namespace roboptim;

typedef boost::mpl::list< ::roboptim::EigenMatrixDense,
  ::roboptim::EigenMatrixSparse> functionTypes_t;

struct F : public DifferentiableFunction
{
  F () : DifferentiableFunction (22, 10, "f_n (x) = n * x")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < outputSize (); ++i)
      for (size_type j = 0; j < 3; ++j)
	{
	  res[i] += (value_type)i * argument[19 + j];
	}
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
		      size_type functionId) const
  {
    grad.setZero ();
    for (size_type j = 0; j < 3; ++j)
      {
	grad[19 + j] += (value_type)functionId;
      }
  }
};

struct G : public DifferentiableFunction
{
  G () : DifferentiableFunction (3, 1, "f_n (x) = sum(x)")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < inputSize (); ++i)
      {
	res[0] += argument[i];
      }
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
		      size_type functionId) const
  {
    grad.setZero ();
    for (size_type j = 0; j < 3; ++j)
      grad[j] = functionId * 0 + 1;
  }
};

struct H : public DifferentiableFunction
{
  H () : DifferentiableFunction (45, 3, "f_n (x) = sum(x)")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < inputSize (); ++i)
      {
	res[i % 3] += argument[i];
      }
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
		      size_type functionId) const
  {
    grad.setZero ();
    for (size_type j = 0; j < 45; ++j)
      grad[j] = (functionId == (j % 3));
  }
};

boost::shared_ptr<boost::test_tools::output_test_stream> output;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (manifold_map_test_0)
{/*
  output = retrievePattern("filter-manifold-map");

  boost::shared_ptr<F> f (new F());

  pgs::RealSpace pos(3);pos.name() = "position";
  pgs::SO3<pgs::ExpMapMatrix> ori; ori.name() = "orientation";
  const pgs::CartesianProduct freeFlyer(pos, ori);
  pgs::RealSpace joints(10); joints.name() = "joints";
  const pgs::CartesianProduct robot(freeFlyer, joints);

  const pgs::CartesianProduct cartProd(joints, ori);
  const pgs::CartesianProduct myFuncManifold(cartProd, pos);

  boost::shared_ptr<DescriptiveWrapper<DifferentiableFunction>>
    descWrapPtr(new DescriptiveWrapper<DifferentiableFunction>(f, myFuncManifold));

  InstanceWrapper<DifferentiableFunction> instWrap(descWrapPtr, robot, myFuncManifold);

  InstanceWrapper<DifferentiableFunction>::argument_t input = Eigen::VectorXd::Zero(22);
  InstanceWrapper<DifferentiableFunction>::result_t result = Eigen::VectorXd::Zero(10);
  InstanceWrapper<DifferentiableFunction>::gradient_t gradient = Eigen::VectorXd::Zero(22);
  InstanceWrapper<DifferentiableFunction>::jacobian_t jacobian = Eigen::MatrixXd::Zero(10, 22);

  for(int i = 0; i < 3; ++i)
    {
      input(i) = 1 + i;
    }

  (*output) << instWrap << "\n";
  std::cout << instWrap << std::endl;
  std::cout << "input: " << input.transpose() << std::endl;
  instWrap(result, input);
  std::cout << "result: " << result.transpose() << std::endl << std::endl;

  for (int i = 0; i < f->outputSize(); ++i)
    {
      instWrap.gradient(gradient, input, i);
    }

  instWrap.jacobian(jacobian, input);
  std::cout << "jacobian: " << std::endl << jacobian << std::endl;

  (*output) << (*descWrapPtr);

  BOOST_CHECK (output->match_pattern());

}

BOOST_AUTO_TEST_CASE (manifold_map_test_1)
{
  output = retrievePattern("filter-manifold-map-1");

  boost::shared_ptr<G> g (new G());

  std::vector<const pgs::RealSpace*> reals;
  pgs::CartesianProduct problemManifold;
  const pgs::RealSpace descriptiveManifold(3);

  boost::shared_ptr<DescriptiveWrapper<DifferentiableFunction>>
    descWrapPtr(new DescriptiveWrapper<DifferentiableFunction>(g, descriptiveManifold));

  size_t posNumber = 15;

  for (size_t i = 0; i < posNumber; ++i)
    {
      pgs::RealSpace* newR = new pgs::RealSpace(3);
      newR->name() = "position (" + std::to_string(i) + ")";
      reals.push_back(newR);
      problemManifold.multiply(*reals.back());
    }

  InstanceWrapper<DifferentiableFunction>::argument_t input = Eigen::VectorXd::Zero(3 * static_cast<long>(posNumber));
  InstanceWrapper<DifferentiableFunction>::result_t result = Eigen::VectorXd::Zero(1);
  InstanceWrapper<DifferentiableFunction>::gradient_t gradient = Eigen::VectorXd::Zero(3 * static_cast<long>(posNumber));
  InstanceWrapper<DifferentiableFunction>::jacobian_t jacobian = Eigen::MatrixXd::Zero(1, 3 * static_cast<long>(posNumber));

  for (int i = 0; i < input.size(); ++i)
    {
      input(i) = i;
    }

  for (size_t i = 0; i < posNumber; ++i)
    {
      InstanceWrapper<DifferentiableFunction> instWrap(descWrapPtr, problemManifold, *reals[i]);

      instWrap(result, input);

      BOOST_CHECK(result(0) == (3 * (3 * i + 1)));

      instWrap.jacobian(jacobian, input);
      std::cout << jacobian << std::endl;
      (*output) << jacobian << std::endl;
    }

  for (size_t i = 0; i < posNumber; ++i)
    {
      delete reals[i];
    }

  }

BOOST_AUTO_TEST_CASE (manifold_map_test_2)
{
  output = retrievePattern("filter-manifold-map-2");

  boost::shared_ptr<F> f (new F());
  boost::shared_ptr<G> g (new G());

  const pgs::RealSpace manifold(2);

  try
  {
    boost::shared_ptr<DescriptiveWrapper<DifferentiableFunction>>
    descWrapPtr(new DescriptiveWrapper<DifferentiableFunction>(g, manifold));
  }
  catch (std::runtime_error& e)
  {
    (*output) << "std::runtime_error: " << e.what() << "\n";
  }

  const pgs::SO3<pgs::ExpMapMatrix> manifold2;
  try
  {
    boost::shared_ptr<DescriptiveWrapper<DifferentiableFunction>>
    descWrapPtr(new DescriptiveWrapper<DifferentiableFunction>(f, manifold2));
  }
  catch (std::runtime_error& e)
  {
    (*output) << "std::runtime_error: " << e.what() << "\n";
  }

  const pgs::CartesianProduct manifold3(manifold, manifold2);

  try
  {
    boost::shared_ptr<DescriptiveWrapper<DifferentiableFunction>>
    descWrapPtr(new DescriptiveWrapper<DifferentiableFunction>(g, manifold3));
  }
  catch (std::runtime_error& e)
  {
    (*output) << "std::runtime_error: " << e.what() << "\n";
  }

  BOOST_CHECK (output->match_pattern());

}

BOOST_AUTO_TEST_CASE (manifold_map_test_3)
{
  output = retrievePattern("filter-manifold-map-3");

  boost::shared_ptr<H> h (new H());

  std::vector<const pgs::Manifold*> reals;
  std::vector<std::pair<long, long>> restrictions;
  pgs::CartesianProduct problemManifold;
  pgs::CartesianProduct descriptiveManifold;

  const size_t posNumber = 15;

  for (size_t i = 0; i < posNumber; ++i)
    {
      pgs::RealSpace* newR = new pgs::RealSpace(6);
      newR->name() = "position (" + std::to_string(i) + ")";
      reals.push_back(newR);
      problemManifold.multiply(*reals.back());
      descriptiveManifold.multiply(*(new pgs::RealSpace(3)));
    }

  restrictions.push_back(std::make_pair(3l, 3l));

  boost::shared_ptr<DescriptiveWrapper<DifferentiableFunction>>
    descWrapPtr(new DescriptiveWrapper<DifferentiableFunction>(h, descriptiveManifold));

  InstanceWrapper<DifferentiableFunction>::argument_t input = Eigen::VectorXd::Zero(6 * static_cast<long>(posNumber));
  InstanceWrapper<DifferentiableFunction>::result_t result = Eigen::VectorXd::Zero(3);
  InstanceWrapper<DifferentiableFunction>::gradient_t gradient = Eigen::VectorXd::Zero(6 * static_cast<long>(posNumber));
  InstanceWrapper<DifferentiableFunction>::jacobian_t jacobian = Eigen::MatrixXd::Zero(3, 6 * static_cast<long>(posNumber));

  InstanceWrapper<DifferentiableFunction> instWrap(descWrapPtr, problemManifold, problemManifold, reals, restrictions);

 for (int i = 0; i < input.size(); ++i)
    {
      input(i) = (i % 6) > 2;
    }

 instWrap(result, input);

 (*output) << instWrap;

 std::cout << "result: " << result << std::endl;

 BOOST_CHECK (result(0) == 15);
 BOOST_CHECK (result(1) == 15);
 BOOST_CHECK (result(2) == 15);

 instWrap.jacobian(jacobian, input);

 Eigen::MatrixXd jac = jacobian;

 (*output) << "\n";

 for (int i = 0; i < jac.rows(); ++i)
   {
     std::cout << jac.row(i) << std::endl;
     (*output) << jac.row(i) << std::endl;
   }

 BOOST_CHECK (output->match_pattern());
}

BOOST_AUTO_TEST_CASE (manifold_map_test_4)
{
  boost::shared_ptr<F> f (new F());

  pgs::RealSpace pos(3);pos.name() = "position";
  pgs::SO3<pgs::ExpMapMatrix> ori; ori.name() = "orientation";
  const pgs::CartesianProduct freeFlyer(pos, ori);
  pgs::RealSpace joints(10); joints.name() = "joints";
  const pgs::CartesianProduct robot(freeFlyer, joints);

  const pgs::S2 s2;
  const pgs::CartesianProduct cartProd(joints, ori);
  const pgs::CartesianProduct myFuncManifold(cartProd, s2);
  const pgs::CartesianProduct mySubManifold(cartProd, pos);

  boost::shared_ptr<DescriptiveWrapper<DifferentiableFunction>>
    descWrapPtr(new DescriptiveWrapper<DifferentiableFunction>(f, myFuncManifold));

  std::vector<const pgs::Manifold*> restrictedManifolds;
  restrictedManifolds.push_back(&pos);
  std::vector<std::pair<long, long>> restrictions;
  restrictions.push_back(std::make_pair(1l, 1l));

  bool errorThrown = false;

  try
    {
      InstanceWrapper<DifferentiableFunction> instWrap(descWrapPtr, robot, mySubManifold);
    }
  catch (std::runtime_error& e)
    {
      errorThrown = true;
    }

  BOOST_CHECK(errorThrown);
  errorThrown = false;

    try
    {
      InstanceWrapper<DifferentiableFunction> instWrap(descWrapPtr, robot, mySubManifold, restrictedManifolds, restrictions);
    }
  catch (std::runtime_error& e)
    {
      errorThrown = true;
    }

  BOOST_CHECK(errorThrown);
*/}


BOOST_AUTO_TEST_SUITE_END ()
