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

typedef boost::mpl::list< ::roboptim::EigenMatrixDense,
			  ::roboptim::EigenMatrixSparse> functionTypes_t;

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

template <>
inline void
F<roboptim::EigenMatrixSparse>::impl_gradient(gradient_ref grad, const_argument_ref, size_type functionId) const
{
  grad.setZero ();
  for (size_type j = 0; j < 3; ++j)
    {
      grad.insert (19 + j) += (value_type)functionId;
    }
}

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

template <>
inline void
G<roboptim::EigenMatrixSparse>::impl_gradient(gradient_ref grad, const_argument_ref, size_type functionId) const
{
  grad.setZero ();
  for (size_type j = 0; j < 3; ++j)
    {
      grad.insert (j) = functionId * 0 + 1;
    }
}

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

template <>
inline void
H<roboptim::EigenMatrixSparse>::impl_gradient(gradient_ref grad, const_argument_ref, size_type functionId) const
{
  grad.setZero ();
  for (size_type j = 0; j < 45; ++j)
    {
      grad.insert (j) = (functionId == (j % 3));
    }
}

template<class T>
struct I : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  I () : roboptim::GenericDifferentiableFunction<T> (15, 15, "f_n (x) = sum(x)")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < this->inputSize (); ++i)
      {
	res[i] = argument[i];
      }
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
		      size_type functionId) const
  {
    grad.setZero ();
    for (size_type j = 0; j < this->inputSize (); ++j)
      grad[j] = (functionId == j?1:0);
  }
};

template <>
inline void
I<roboptim::EigenMatrixSparse>::impl_gradient(gradient_ref grad, const_argument_ref, size_type functionId) const
{
  grad.setZero ();
  for (size_type j = 0; j < this->inputSize (); ++j)
    {
      grad.insert (j) = (functionId == j?1:0);
    }
}

inline roboptim::GenericFunctionTraits<roboptim::EigenMatrixDense>::matrix_t
to_dense (roboptim::GenericFunctionTraits<roboptim::EigenMatrixDense>::const_matrix_ref m)
{
  return m;
}

inline roboptim::GenericFunctionTraits<roboptim::EigenMatrixDense>::matrix_t
to_dense (roboptim::GenericFunctionTraits<roboptim::EigenMatrixSparse>::const_matrix_ref m)
{
  return roboptim::sparse_to_dense (m);
}

boost::shared_ptr<boost::test_tools::output_test_stream> output;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_0, T, functionTypes_t)
{
  output = retrievePattern("manifold-map");

  typedef F<T> Func;

  DESC_MANIFOLD(FreeFlyerPlus10, REAL_SPACE(10), roboptim::SO3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_FreeFlyerPlus10, Func, FreeFlyerPlus10);
  typedef roboptim::FunctionOnManifold<typename Func::parent_t> Instance_F_On_FreeFlyerPlus10;

  mnf::RealSpace pos(3);pos.name() = "position";
  mnf::SO3<mnf::ExpMapMatrix> ori; ori.name() = "orientation";
  const mnf::CartesianProduct freeFlyer(pos, ori);
  mnf::RealSpace joints(10); joints.name() = "joints";
  const mnf::CartesianProduct robot(freeFlyer, joints);

  const mnf::CartesianProduct cartProd(joints, ori);
  const mnf::CartesianProduct myFuncManifold(cartProd, pos);

  F_On_FreeFlyerPlus10 descWrapPtr;

  Instance_F_On_FreeFlyerPlus10 instWrap(descWrapPtr, robot, myFuncManifold);

  typename Instance_F_On_FreeFlyerPlus10::argument_t input = Eigen::VectorXd::Zero(22);
  typename Instance_F_On_FreeFlyerPlus10::result_t result = Eigen::VectorXd::Zero(10);
  typename Instance_F_On_FreeFlyerPlus10::gradient_t gradient(22);
  gradient.setZero();
  typename Instance_F_On_FreeFlyerPlus10::jacobian_t jacobian(10, 22);
  jacobian.setZero();

  Eigen::MatrixXd refJacobian = Eigen::MatrixXd::Zero(10, 22);

  for(int i = 0; i < 3; ++i)
    {
      input(i) = 1 + i;
    }

  (*output) << instWrap << roboptim::iendl;
  instWrap(result, input);

  for (int i = 0; i < 10; ++i)
    {
      instWrap.gradient(gradient, input, i);
    }

  instWrap.jacobian(jacobian, input);
  try
    {
      instWrap.manifold_jacobian(refJacobian, input);
    }
  catch(std::runtime_error e)
    {
      std::cout << e.what() << std::endl;
    }
  (*output) << descWrapPtr << roboptim::iendl;
  std::cout << output->str() << std::endl;
  BOOST_CHECK (output->match_pattern());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_1, T, functionTypes_t)
{
  output = retrievePattern("manifold-map-1");

  typedef G<T> Func;

  DESC_MANIFOLD(Real3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_Real3, Func, Real3);
  typedef roboptim::FunctionOnManifold<typename Func::parent_t> Instance_F_On_Real3;

  std::vector<const mnf::RealSpace*> reals;
  mnf::CartesianProduct problemManifold;
  const mnf::RealSpace descriptiveManifold(3);

  F_On_Real3 descWrapPtr;

  size_t posNumber = 15;

  for (size_t i = 0; i < posNumber; ++i)
    {
      mnf::RealSpace* newR = new mnf::RealSpace(3);
      newR->name() = "position (" + std::to_string(i) + ")";
      reals.push_back(newR);
      problemManifold.multiply(*reals.back());
    }

  typename Instance_F_On_Real3::argument_t input = Eigen::VectorXd::Zero(3 * static_cast<long>(posNumber));
  typename Instance_F_On_Real3::result_t result = Eigen::VectorXd::Zero(1);
  typename Instance_F_On_Real3::gradient_t gradient(3 * static_cast<int> (posNumber));
  gradient.setZero();
  typename Instance_F_On_Real3::jacobian_t jacobian(1, 3 * static_cast<int> (posNumber));
  jacobian.setZero();

  for (int i = 0; i < input.size(); ++i)
    {
      input(i) = i;
    }
  for (size_t i = 0; i < posNumber; ++i)
    {
      Instance_F_On_Real3 instWrap(descWrapPtr, problemManifold, *reals[i]);

      (*output) << "This is the FunctionOnManifold number " << i << ":"
		<< roboptim::incindent;
      (*output) << roboptim::iendl << instWrap;
      instWrap(result, input);
      (*output) << roboptim::iendl << "Result:";
      (*output) << roboptim::iendl << result;

      BOOST_CHECK(result(0) == (3 * (3 * i + 1)));

      instWrap.jacobian(jacobian, input);
      (*output) << roboptim::iendl << "Jacobian:";
      (*output) << roboptim::iendl << to_dense(jacobian) << roboptim::decindent
		<< roboptim::iendl;
    }

  for (size_t i = 0; i < posNumber; ++i)
    {
      delete reals[i];
    }

  std::cout << output->str() << std::endl;
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
      (*output) << "std::runtime_error: " << e.what() << roboptim::iendl;
    }

  DESC_MANIFOLD(Manifold2, roboptim::SO3);
  typedef roboptim::DescriptiveWrapper<Func, Manifold2> F_On_Manifold2;

  try
    {
      new F_On_Manifold2();
    }
  catch (std::runtime_error& e)
    {
      (*output) << "std::runtime_error: " << e.what() << roboptim::iendl;
    }

  DESC_MANIFOLD(Manifold3, REAL_SPACE(2), roboptim::SO3);
  typedef roboptim::DescriptiveWrapper<Gunc, Manifold3> G_On_Manifold3;

  try
    {
      new G_On_Manifold3();
    }
  catch (std::runtime_error& e)
    {
      (*output) << "std::runtime_error: " << e.what() << roboptim::iendl;
    }

  BOOST_CHECK (output->match_pattern());
}

const size_t posNumber = 15;

DEFINE_MANIFOLD(MultipleReal3)
{
  mnf::CartesianProduct* cartesian = new mnf::CartesianProduct();

  for(size_t i = 0; i < posNumber; ++i)
    {
      cartesian->multiply(*(new mnf::RealSpace(3)));
    }

  return cartesian;
}

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_3, T, functionTypes_t)
{
  output = retrievePattern("manifold-map-3");

  typedef H<T> Func;

  std::vector<const mnf::Manifold*> reals;
  std::vector<std::pair<long, long>> restrictions;
  mnf::CartesianProduct problemManifold;
  mnf::CartesianProduct descriptiveManifold;

  for (size_t i = 0; i < posNumber; ++i)
    {
      mnf::RealSpace* newR = new mnf::RealSpace(6);
      newR->name() = "position (" + std::to_string(i) + ")";
      reals.push_back(newR);
      problemManifold.multiply(*reals.back());
      descriptiveManifold.multiply(*(new mnf::RealSpace(3)));
    }

  restrictions.push_back(std::make_pair(3l, 3l));

  DESC_MANIFOLD(MultipleReal3, Manifold_MultipleReal3);
  NAMED_FUNCTION_BINDING(H_On_MultipleReal3, Func, MultipleReal3);
  typedef roboptim::FunctionOnManifold<typename Func::parent_t> Instance_H_On_MultipleReal3;

  H_On_MultipleReal3 descWrapPtr;

  typename Instance_H_On_MultipleReal3::argument_t input = Eigen::VectorXd::Zero(6 * static_cast<long>(posNumber));
  typename Instance_H_On_MultipleReal3::result_t result = Eigen::VectorXd::Zero(3);
  typename Instance_H_On_MultipleReal3::gradient_t gradient(6 * static_cast<int> (posNumber));
  gradient.setZero();
  typename Instance_H_On_MultipleReal3::jacobian_t jacobian(3, 6 * static_cast<int> (posNumber));
  jacobian.setZero();

  Instance_H_On_MultipleReal3 instWrap(descWrapPtr, problemManifold, problemManifold, reals, restrictions);

  for (int i = 0; i < input.size(); ++i)
    {
      input(i) = (i % 6) > 2;
    }

  instWrap(result, input);

  (*output) << instWrap;

  BOOST_CHECK (result(0) == 15);
  BOOST_CHECK (result(1) == 15);
  BOOST_CHECK (result(2) == 15);

  instWrap.jacobian(jacobian, input);

  Eigen::MatrixXd jac = to_dense(jacobian);

  (*output) << roboptim::iendl;

  for (int i = 0; i < jac.rows(); ++i)
    {
      (*output) << jac.row(i) << roboptim::iendl;
    }

  BOOST_CHECK (output->match_pattern());
}

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_test_4, T, functionTypes_t)
{
  typedef F<T> Func;

  mnf::RealSpace pos(3);pos.name() = "position";
  mnf::SO3<mnf::ExpMapMatrix> ori; ori.name() = "orientation";
  const mnf::CartesianProduct freeFlyer(pos, ori);
  mnf::RealSpace joints(10); joints.name() = "joints";
  const mnf::CartesianProduct robot(freeFlyer, joints);

  const mnf::S2 s2;
  const mnf::CartesianProduct cartProd(joints, ori);
  const mnf::CartesianProduct myFuncManifold(cartProd, s2);
  const mnf::CartesianProduct mySubManifold(cartProd, pos);

  DESC_MANIFOLD(FreeFlyerPlus10, REAL_SPACE(10), roboptim::SO3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_FreeFlyerPlus10, Func, FreeFlyerPlus10);
  typedef roboptim::FunctionOnManifold<typename Func::parent_t> Instance_F_On_FreeFlyerPlus10;

  F_On_FreeFlyerPlus10 descWrapPtr;

  std::vector<const mnf::Manifold*> restrictedManifolds;
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

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_split_manifold_into_pieces, T, functionTypes_t)
{
  typedef I<T> Iunc;

  DESC_MANIFOLD(R5X3, REAL_SPACE(5), REAL_SPACE(5), REAL_SPACE(5));
  NAMED_FUNCTION_BINDING(I_On_R5X3, Iunc, R5X3);
  typedef roboptim::FunctionOnManifold<typename Iunc::parent_t> Instance_I_On_R5X3;

  mnf::RealSpace r15(15);
  mnf::CartesianProduct splitR15;
  splitR15.multiply(r15).multiply(r15).multiply(r15);

  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;

  restricted.push_back(&r15);
  restricted.push_back(&r15);
  restricted.push_back(&r15);

  restrictions.push_back(std::make_pair(5l, 5l));
  restrictions.push_back(std::make_pair(10l, 5l));
  restrictions.push_back(std::make_pair(0l, 5l));

  I_On_R5X3 descI;
  Instance_I_On_R5X3 funcI(descI, r15, splitR15, restricted, restrictions);

  typename Instance_I_On_R5X3::argument_t input = Eigen::VectorXd::Zero(15);
  typename Instance_I_On_R5X3::result_t result = Eigen::VectorXd::Zero(15);

  for(int i = 0; i < input.size(); ++i)
    {
      input(i) = i;
    }

  funcI(result, input);

  for(int i = 0; i < 15; ++i)
    {
      BOOST_CHECK(result(i) == ((static_cast<int>(input(i)) + 5) % 15));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_split_fail_not_enough_restrictions, T, functionTypes_t)
{
  typedef I<T> Iunc;

  DESC_MANIFOLD(R5X3, REAL_SPACE(5), REAL_SPACE(5), REAL_SPACE(5));
  NAMED_FUNCTION_BINDING(I_On_R5X3, Iunc, R5X3);
  typedef roboptim::FunctionOnManifold<typename Iunc::parent_t> Instance_I_On_R5X3;

  mnf::RealSpace r15(15);
  mnf::CartesianProduct splitR15;
  splitR15.multiply(r15).multiply(r15).multiply(r15);

  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;

  restricted.push_back(&r15);
  restricted.push_back(&r15);

  restrictions.push_back(std::make_pair(5l, 5l));
  restrictions.push_back(std::make_pair(10l, 5l));

  I_On_R5X3 descI;

  bool errorThrown = false;

  try
    {
      Instance_I_On_R5X3 funcI(descI, r15, splitR15, restricted, restrictions);
    }
  catch (std::runtime_error& e)
    {
      errorThrown = true;
    }

  BOOST_CHECK(errorThrown);
  errorThrown = false;

}

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_map_split_fail_restriction_on_unknown_manifold, T, functionTypes_t)
{
  typedef I<T> Iunc;

  DESC_MANIFOLD(R5X3, REAL_SPACE(5), REAL_SPACE(5), REAL_SPACE(5));
  NAMED_FUNCTION_BINDING(I_On_R5X3, Iunc, R5X3);
  typedef roboptim::FunctionOnManifold<typename Iunc::parent_t> Instance_I_On_R5X3;

  mnf::RealSpace r15(15);
  mnf::CartesianProduct splitR15;
  mnf::SO3<mnf::ExpMapMatrix> so3;
  splitR15.multiply(r15).multiply(r15).multiply(r15);

  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;

  restricted.push_back(&r15);
  restricted.push_back(&r15);
  restricted.push_back(&r15);
  restricted.push_back(&so3);

  restrictions.push_back(std::make_pair(5l, 5l));
  restrictions.push_back(std::make_pair(10l, 5l));
  restrictions.push_back(std::make_pair(0l, 5l));
  restrictions.push_back(std::make_pair(3l, 3l));

  I_On_R5X3 descI;

  bool errorThrown = false;

  try
    {
      Instance_I_On_R5X3 funcI(descI, r15, splitR15, restricted, restrictions);
    }
  catch (std::runtime_error& e)
    {
      errorThrown = true;
    }

  BOOST_CHECK(errorThrown);
  errorThrown = false;

}

BOOST_AUTO_TEST_SUITE_END ()
