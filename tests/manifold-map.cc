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

#include <roboptim/core/manifold-map/decorator/manifold-map.hh>
#include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>
#include <roboptim/core/manifold-map/decorator/problem-factory.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>

//using namespace roboptim;

typedef boost::mpl::list< ::roboptim::EigenMatrixDense/*,
			  ::roboptim::EigenMatrixSparse*/> functionTypes_t;

template<class T>
struct F : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  F () : roboptim::GenericDifferentiableFunction<T> (22, 10, "f_n (x) = n * x")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < this->outputSize (); ++i)
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

template<class T>
struct G : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  G () : roboptim::GenericDifferentiableFunction<T> (3, 1, "f_n (x) = sum(x)")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < this->inputSize (); ++i)
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

template<class T>
struct H : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  H () : roboptim::GenericDifferentiableFunction<T> (45, 3, "f_n (x) = sum(x)")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < this->inputSize (); ++i)
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

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_0, T, functionTypes_t)
{
  output = retrievePattern("manifold-map");

  typedef F<T> Func;

  DESC_MANIFOLD(FreeFlyerPlus10, REAL_SPACE(10), roboptim::SO3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_FreeFlyerPlus10, Func, FreeFlyerPlus10);

  pgs::RealSpace pos(3);pos.name() = "position";
  pgs::SO3<pgs::ExpMapMatrix> ori; ori.name() = "orientation";
  const pgs::CartesianProduct freeFlyer(pos, ori);
  pgs::RealSpace joints(10); joints.name() = "joints";
  const pgs::CartesianProduct robot(freeFlyer, joints);

  const pgs::CartesianProduct cartProd(joints, ori);
  const pgs::CartesianProduct myFuncManifold(cartProd, pos);

  F_On_FreeFlyerPlus10 descWrapPtr;

  Instance_F_On_FreeFlyerPlus10 instWrap(descWrapPtr, robot, myFuncManifold);

  typename Instance_F_On_FreeFlyerPlus10::argument_t input = Eigen::VectorXd::Zero(22);
  typename Instance_F_On_FreeFlyerPlus10::result_t result = Eigen::VectorXd::Zero(10);
  typename Instance_F_On_FreeFlyerPlus10::gradient_t gradient = Eigen::VectorXd::Zero(22);
  typename Instance_F_On_FreeFlyerPlus10::jacobian_t jacobian = Eigen::MatrixXd::Zero(10, 22);

  Eigen::MatrixXd refJacobian = Eigen::MatrixXd::Zero(10, 22);

  for(int i = 0; i < 3; ++i)
    {
      input(i) = 1 + i;
    }

  (*output) << instWrap << "\n";
  std::cout << instWrap << std::endl;
  std::cout << "input: " << input.transpose() << std::endl;
  instWrap(result, input);
  std::cout << "result: " << result.transpose() << std::endl << std::endl;

  for (int i = 0; i < 10; ++i)
    {
      std::cout << "i: " << i << std::endl;
      instWrap.gradient(gradient, input, i);
    }

  instWrap.jacobian(jacobian, input);
  std::cout << "jacobian: " << std::endl << jacobian << std::endl;
  instWrap.manifold_jacobian(refJacobian, input);

  (*output) << descWrapPtr;

  BOOST_CHECK (output->match_pattern());

 }

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_1, T, functionTypes_t)
{
  output = retrievePattern("manifold-map-1");

  typedef G<T> Func;

  DESC_MANIFOLD(Real3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_Real3, Func, Real3);

  std::vector<const pgs::RealSpace*> reals;
  pgs::CartesianProduct problemManifold;
  const pgs::RealSpace descriptiveManifold(3);

 F_On_Real3 descWrapPtr;

  size_t posNumber = 15;

  for (size_t i = 0; i < posNumber; ++i)
    {
      pgs::RealSpace* newR = new pgs::RealSpace(3);
      newR->name() = "position (" + std::to_string(i) + ")";
      reals.push_back(newR);
      problemManifold.multiply(*reals.back());
    }

  typename Instance_F_On_Real3::argument_t input = Eigen::VectorXd::Zero(3 * static_cast<long>(posNumber));
  typename Instance_F_On_Real3::result_t result = Eigen::VectorXd::Zero(1);
  typename Instance_F_On_Real3::gradient_t gradient = Eigen::VectorXd::Zero(3 * static_cast<long>(posNumber));
  typename Instance_F_On_Real3::jacobian_t jacobian = Eigen::MatrixXd::Zero(1, 3 * static_cast<long>(posNumber));

  for (int i = 0; i < input.size(); ++i)
    {
      input(i) = i;
    }

  for (size_t i = 0; i < posNumber; ++i)
    {
      Instance_F_On_Real3 instWrap(descWrapPtr, problemManifold, *reals[i]);

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

    BOOST_CHECK (output->match_pattern());
  }

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_2, T, functionTypes_t)
{
  output = retrievePattern("manifold-map-2");

  typedef F<T> Func;
  typedef G<T> Gunc;

  DESC_MANIFOLD(Real2, REAL_SPACE(2));
  typedef roboptim::DescriptiveWrapper<Gunc, Real2> G_On_Real2;

  try
  {
    new G_On_Real2();
  }
  catch (std::runtime_error& e)
  {
    (*output) << "std::runtime_error: " << e.what() << "\n";
  }

  DESC_MANIFOLD(Manifold2, roboptim::SO3);
  typedef roboptim::DescriptiveWrapper<Func, Manifold2> F_On_Manifold2;

  try
  {
    new F_On_Manifold2();
  }
  catch (std::runtime_error& e)
  {
    (*output) << "std::runtime_error: " << e.what() << "\n";
  }

  DESC_MANIFOLD(Manifold3, REAL_SPACE(2), roboptim::SO3);
  typedef roboptim::DescriptiveWrapper<Func, Manifold3> F_On_Manifold3;

  try
  {
    new F_On_Manifold3();
  }
  catch (std::runtime_error& e)
  {
    (*output) << "std::runtime_error: " << e.what() << "\n";
  }

  BOOST_CHECK (output->match_pattern());

}

const size_t posNumber = 15;

DEFINE_MANIFOLD(MultipleReal3)
{
  pgs::CartesianProduct* cartesian = new pgs::CartesianProduct();

  for(size_t i = 0; i < posNumber; ++i)
    {
      cartesian->multiply(*(new pgs::RealSpace(3)));
    }

  return cartesian;
}

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_3, T, functionTypes_t)
{
  output = retrievePattern("manifold-map-3");

  typedef H<T> Func;

  std::vector<const pgs::Manifold*> reals;
  std::vector<std::pair<long, long>> restrictions;
  pgs::CartesianProduct problemManifold;
  pgs::CartesianProduct descriptiveManifold;

  for (size_t i = 0; i < posNumber; ++i)
    {
      pgs::RealSpace* newR = new pgs::RealSpace(6);
      newR->name() = "position (" + std::to_string(i) + ")";
      reals.push_back(newR);
      problemManifold.multiply(*reals.back());
      descriptiveManifold.multiply(*(new pgs::RealSpace(3)));
    }

  restrictions.push_back(std::make_pair(3l, 3l));

  DESC_MANIFOLD(MultipleReal3, Manifold_MultipleReal3);
  NAMED_FUNCTION_BINDING(H_On_MultipleReal3, Func, MultipleReal3);

  H_On_MultipleReal3 descWrapPtr;

  typename Instance_H_On_MultipleReal3::argument_t input = Eigen::VectorXd::Zero(6 * static_cast<long>(posNumber));
  typename Instance_H_On_MultipleReal3::result_t result = Eigen::VectorXd::Zero(3);
  typename Instance_H_On_MultipleReal3::gradient_t gradient = Eigen::VectorXd::Zero(6 * static_cast<long>(posNumber));
  typename Instance_H_On_MultipleReal3::jacobian_t jacobian = Eigen::MatrixXd::Zero(3, 6 * static_cast<long>(posNumber));

  Instance_H_On_MultipleReal3 instWrap(descWrapPtr, problemManifold, problemManifold, reals, restrictions);

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

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_4, T, functionTypes_t)
{
  typedef F<T> Func;

  pgs::RealSpace pos(3);pos.name() = "position";
  pgs::SO3<pgs::ExpMapMatrix> ori; ori.name() = "orientation";
  const pgs::CartesianProduct freeFlyer(pos, ori);
  pgs::RealSpace joints(10); joints.name() = "joints";
  const pgs::CartesianProduct robot(freeFlyer, joints);

  const pgs::S2 s2;
  const pgs::CartesianProduct cartProd(joints, ori);
  const pgs::CartesianProduct myFuncManifold(cartProd, s2);
  const pgs::CartesianProduct mySubManifold(cartProd, pos);

  DESC_MANIFOLD(FreeFlyerPlus10, REAL_SPACE(10), roboptim::SO3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_FreeFlyerPlus10, Func, FreeFlyerPlus10);

  F_On_FreeFlyerPlus10 descWrapPtr;

  std::vector<const pgs::Manifold*> restrictedManifolds;
  restrictedManifolds.push_back(&pos);
  std::vector<std::pair<long, long>> restrictions;
  restrictions.push_back(std::make_pair(1l, 1l));

  bool errorThrown = false;

  try
    {
      Instance_F_On_FreeFlyerPlus10 instWrap(descWrapPtr, robot, myFuncManifold);
    }
  catch (std::runtime_error& e)
    {
      errorThrown = true;
    }

  BOOST_CHECK(errorThrown);
  errorThrown = false;

    try
    {
      Instance_F_On_FreeFlyerPlus10 instWrap(descWrapPtr, robot, mySubManifold, restrictedManifolds, restrictions);
    }
  catch (std::runtime_error& e)
    {
      errorThrown = true;
    }

  BOOST_CHECK(errorThrown);
}

BOOST_AUTO_TEST_SUITE_END ()
