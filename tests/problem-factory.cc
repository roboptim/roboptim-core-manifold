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
#include <roboptim/core/manifold-map/decorator/manifold-problem-factory.hh>

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

  F () : roboptim::GenericDifferentiableFunction<T> (9, 1, "F(x) = sum(x)")
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
		      size_type) const
  {
    grad.setOnes ();
  }

  virtual ~F()
  {
    std::cout << "END OF F" << std::endl;
  }
};

template<>
inline void
F<roboptim::EigenMatrixSparse>::impl_gradient (gradient_ref grad, const_argument_ref,
					       size_type) const
{
  for (size_type i =0; i < this->inputSize(); ++i)
    {
      grad.coeffRef(i) = 1;
    }
}

template<class T>
struct G : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  G () : roboptim::GenericDifferentiableFunction<T> (3, 1, "G(x) = sum(x)")
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
		      size_type) const
  {
    grad.setOnes ();
  }
};

template<>
inline void
G<roboptim::EigenMatrixSparse>::impl_gradient (gradient_ref grad, const_argument_ref,
					       size_type) const
{
  for (size_type i =0; i < this->inputSize(); ++i)
    {
      grad.coeffRef(i) = 1;
    }
}

template<class T>
struct H : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  H () : roboptim::GenericDifferentiableFunction<T> (10, 4, "H(x) = sum(x)")
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
		      size_type) const
  {
    grad.setOnes ();
  }
};

template<>
inline void
H<roboptim::EigenMatrixSparse>::impl_gradient (gradient_ref grad, const_argument_ref,
					       size_type) const
{
  for (size_type i =0; i < this->inputSize(); ++i)
    {
      grad.coeffRef(i) = 1;
    }
}

template<class T>
struct I : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  I () : roboptim::GenericDifferentiableFunction<T> (22, 1, "I(x) = sum(x)")
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
		      size_type) const
  {
    grad.setOnes ();
  }
};

template<>
inline void
I<roboptim::EigenMatrixSparse>::impl_gradient (gradient_ref grad, const_argument_ref,
					       size_type) const
{
  for (size_type i =0; i < this->inputSize(); ++i)
    {
      grad.coeffRef(i) = 1;
    }
}


boost::shared_ptr<boost::test_tools::output_test_stream> output;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_factory_test, T, functionTypes_t)
{
  typedef F<T> Func;
  typedef G<T> Gunc;
  typedef H<T> Hunc;
  typedef I<T> Iunc;

  roboptim::ManifoldProblemFactory<T> factory;

  ROBOPTIM_DESC_MANIFOLD(R3, ROBOPTIM_REAL_SPACE(3));
  ROBOPTIM_NAMED_FUNCTION_BINDING(G_On_R3, Gunc, R3);

  ROBOPTIM_DESC_MANIFOLD(SO3, roboptim::SO3);
  ROBOPTIM_NAMED_FUNCTION_BINDING(F_On_SO3, Func, SO3);

  ROBOPTIM_DESC_MANIFOLD(R10, ROBOPTIM_REAL_SPACE(10));
  ROBOPTIM_NAMED_FUNCTION_BINDING(H_On_R10, Hunc, R10);

  ROBOPTIM_DESC_MANIFOLD(R3XSO3XR10, ROBOPTIM_REAL_SPACE(3), roboptim::SO3, ROBOPTIM_REAL_SPACE(10));
  ROBOPTIM_NAMED_FUNCTION_BINDING(I_On_R3XSO3XR10, Iunc, R3XSO3XR10);

  mnf::RealSpace pos(3);
  mnf::SO3<mnf::ExpMapMatrix> ori;
  mnf::RealSpace joints(10);
  mnf::CartesianProduct prod;
  prod.multiply(pos).multiply(ori).multiply(joints);
  mnf::RealSpace r42(42);
  mnf::RealSpace r39(39);

  mnf::CartesianProduct prod2;
  prod2.multiply(r42).multiply(ori).multiply(r42);
  mnf::CartesianProduct prod3;
  prod3.multiply(r39).multiply(ori).multiply(r39);


  F_On_SO3 cnstr1;
  G_On_R3 objDesc;
  H_On_R10 cnstr2;
  I_On_R3XSO3XR10 cnstr3;

  std::vector<const mnf::Manifold*> restricted;
  std::vector<std::pair<long, long>> restrictions;

  restricted.push_back(&r39);
  restrictions.push_back(std::make_pair(14, 3));
  restricted.push_back(&r39);
  restrictions.push_back(std::make_pair(27, 10));

  factory.addConstraint(cnstr1, ori);
  factory.addConstraint(cnstr2, joints);
  factory.addConstraint(cnstr2, joints);
  factory.addConstraint(objDesc, pos);
  factory.addConstraint(cnstr3, prod3, restricted, restrictions);

  {
    typename Hunc::intervals_t bounds;

    for(int i = 0; i < Hunc().outputSize(); ++i)
      {
	bounds.push_back(roboptim::Function::makeLowerInterval(25));
      }

    factory.addConstraint(cnstr2, joints).setBounds(bounds);
  }

  restricted.clear();
  restricted.push_back(&r42);
  restricted.push_back(&r42);

  factory.setObjective(cnstr3, prod2, restricted, restrictions);

  roboptim::ProblemOnManifold<T>* manifoldProblem = factory.getProblem();

  BOOST_CHECK(manifoldProblem->getManifold().representationDim() == 22 + 42 + 39);

  delete manifoldProblem;
}

BOOST_AUTO_TEST_CASE_TEMPLATE (manifold_factory_no_objective_test, T, functionTypes_t)
{
  typedef F<T> Func;
  typedef G<T> Gunc;
  typedef H<T> Hunc;

  roboptim::ManifoldProblemFactory<T> factory;

  ROBOPTIM_DESC_MANIFOLD(R3, ROBOPTIM_REAL_SPACE(3));
  ROBOPTIM_NAMED_FUNCTION_BINDING(G_On_R3, Gunc, R3);

  ROBOPTIM_DESC_MANIFOLD(SO3, roboptim::SO3);
  ROBOPTIM_NAMED_FUNCTION_BINDING(F_On_SO3, Func, SO3);

  ROBOPTIM_DESC_MANIFOLD(R10, ROBOPTIM_REAL_SPACE(10));
  ROBOPTIM_NAMED_FUNCTION_BINDING(H_On_R10, Hunc, R10);

  mnf::RealSpace pos(3);
  mnf::SO3<mnf::ExpMapMatrix> ori;
  mnf::RealSpace joints(10);

  F_On_SO3 cnstr1;
  G_On_R3 objDesc;
  H_On_R10 cnstr2;

  factory.addConstraint(cnstr1, ori);
  factory.addConstraint(cnstr2, joints);

  {
    std::vector<double> scales;

    for(int i = 0; i < Hunc().outputSize(); ++i)
      {
	scales.push_back(1.);
      }

    factory.addConstraint(cnstr2, joints).setScaling(scales);
  }

  roboptim::ProblemOnManifold<T>* manifoldProblem = factory.getProblem();

  BOOST_CHECK(manifoldProblem->getManifold().representationDim() == 19);

  delete manifoldProblem;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(manifold_problem_factory_no_constraints, T, functionTypes_t)
{
  roboptim::ManifoldProblemFactory<T> factory;

  roboptim::ProblemOnManifold<T>* manifoldProblem;
  BOOST_CHECK_THROW(manifoldProblem = factory.getProblem(); delete manifoldProblem, std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END ()
