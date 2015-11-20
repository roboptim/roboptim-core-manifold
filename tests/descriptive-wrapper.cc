#include "shared-tests/fixture.hh"

#include <roboptim/core/differentiable-function.hh>

#include <roboptim/core/manifold-map/decorator/descriptive-wrapper.hh>
#include <roboptim/core/manifold-map/decorator/manifold-map.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>

typedef boost::mpl::list< ::roboptim::EigenMatrixDense,
			  ::roboptim::EigenMatrixSparse> functionTypes_t;

template <class T>
struct F : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_(
      roboptim::GenericDifferentiableFunction<T>);

  F() : roboptim::GenericDifferentiableFunction<T>(7, 4, "dummy function") {}

  void impl_compute(result_ref, const_argument_ref) const {}

  void impl_gradient(gradient_ref, const_argument_ref, size_type) const {}
};

boost::shared_ptr<boost::test_tools::output_test_stream> output;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (descWrapTest, T, functionTypes_t)
{
  output = retrievePattern("descriptive-wrapper");

  ROBOPTIM_DESC_MANIFOLD(Real7, ROBOPTIM_REAL_SPACE(7));
  ROBOPTIM_DESC_MANIFOLD(Real8, ROBOPTIM_REAL_SPACE(8));
  typedef roboptim::DescriptiveWrapper<F<T>, Real7> FonR7;
  typedef roboptim::DescriptiveWrapper<F<T>, Real8> FonR8;

  {
    (*output) << "DescriptiveWrapper (U* f, const mnf::Manifold& m)"
              << std::endl;
    mnf::RealSpace R7(7);
    mnf::RealSpace R8(8);
    const mnf::Manifold& mR7 = R7; 
    const mnf::Manifold& mR8 = R8; 
    const F<T>* f = new F<T>();
    FonR7 b(f, mR7);
    (*output) << b << std::endl;
    BOOST_CHECK_THROW (FonR7 c(f, mR8), std::runtime_error);   ;
  }
  {
    (*output) << "DescriptiveWrapper (const U* f, const mnf::Manifold& m)"
              << std::endl;
    mnf::RealSpace R7(7);
    mnf::RealSpace R8(8);
    const mnf::Manifold& mR7 = R7; 
    const mnf::Manifold& mR8 = R8; 
    F<T>* f = new F<T>();
    FonR7 b(f, mR7);
    (*output) << b << std::endl;
    BOOST_CHECK_THROW (FonR7 c(f, mR8), std::runtime_error);   ;
  }
  {
    (*output) << "DescriptiveWrapper (std::shared_ptr<const U> f, "
                 "std::shared_ptr<const mnf::Manifold> m)" << std::endl;
    std::shared_ptr<const mnf::Manifold> R7(new mnf::RealSpace(7));
    std::shared_ptr<const mnf::Manifold> R8(new mnf::RealSpace(8));
    std::shared_ptr<const F<T>> f(new F<T>());
    FonR7 b(f, R7);
    (*output) << b << std::endl;
    BOOST_CHECK_THROW (FonR7 b(f, R8), std::runtime_error);   ;
  }
  {
    (*output) << "DescriptiveWrapper (boost::shared_ptr<const U> f, "
                 "boost::shared_ptr<const mnf::Manifold> m)" << std::endl;
    boost::shared_ptr<const mnf::Manifold> R7(new mnf::RealSpace(7));
    boost::shared_ptr<const mnf::Manifold> R8(new mnf::RealSpace(8));
    boost::shared_ptr<const F<T>> f(new F<T>());
    FonR7 b(f, R7);
    (*output) << b << std::endl;
    BOOST_CHECK_THROW (FonR7 b(f, R8), std::runtime_error);   ;
  }
  {
    (*output) << "explicit DescriptiveWrapper (Types ... args)" << std::endl;
    FonR7 a;
    (*output) << a << std::endl;
  }
  {
    (*output) << "makeUNCHECKEDDescriptiveWrapper(const U* fct, const "
                 "mnf::Manifold& manifold)" << std::endl;
    mnf::RealSpace R7(7);
    mnf::RealSpace R8(8);
    F<T>* f = new F<T>();
    auto b = FonR7::makeUNCHECKEDDescriptiveWrapper(f, R7);
    (*output) << *b << std::endl;
    auto c = FonR8::makeUNCHECKEDDescriptiveWrapper(f, R8);
    (*output) << *c << std::endl;
  }
  {
    (*output) << "makeUNCHECKEDDescriptiveWrapper(std::shared_ptr<const U> "
                 "fct, std::shared_ptr<const mnf::Manifold> manifold)"
              << std::endl;
    boost::shared_ptr<mnf::Manifold> R7(new mnf::RealSpace(7));
    boost::shared_ptr<mnf::Manifold> R8(new mnf::RealSpace(8));
    boost::shared_ptr<F<T>> f(new F<T>());
    auto b = FonR7::makeUNCHECKEDDescriptiveWrapper(f, R7);
    (*output) << *b << std::endl;
    auto c = FonR8::makeUNCHECKEDDescriptiveWrapper(f, R8);
    (*output) << *c << std::endl;
  }
  {
    (*output) << "makeUNCHECKEDDescriptiveWrapper( boost::shared_ptr<const U> "
                 "fct, boost::shared_ptr<const mnf::Manifold> manifold)"
              << std::endl;
    std::shared_ptr<mnf::Manifold> R7(new mnf::RealSpace(7));
    std::shared_ptr<mnf::Manifold> R8(new mnf::RealSpace(8));
    std::shared_ptr<F<T>> f(new F<T>());
    auto b = FonR7::makeUNCHECKEDDescriptiveWrapper(f, R7);
    (*output) << *b << std::endl;
    auto c = FonR8::makeUNCHECKEDDescriptiveWrapper(f, R8);
    (*output) << *c << std::endl;
  }
  //std::cout << output->str() << std::endl;
  BOOST_CHECK (output->match_pattern());
}

BOOST_AUTO_TEST_SUITE_END ()
