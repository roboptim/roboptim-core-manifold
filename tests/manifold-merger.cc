#include <iostream>

#include <boost/test/unit_test.hpp>

#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>

#include <roboptim/core/manifold-map/decorator/manifold-merger.hh>

using namespace mnf;
using namespace roboptim;

BOOST_AUTO_TEST_CASE(simpleMerge)
{
  ManifoldMerger mMerger;
  RealSpace R1(1);
  RealSpace R2(2);
  RealSpace R3(3);
  mMerger.addManifold(R1);
  mMerger.addManifold(R2);
  mMerger.addManifold(R3);
  CartesianProduct Res{&R1, &R2, &R3};
  BOOST_CHECK(mMerger.getManifold()->isSameTopology(Res));
  for (size_t i = 0; i < 3; ++i)
    BOOST_CHECK(mMerger.getManifold()->operator()(i).getInstanceId() ==
                Res(i).getInstanceId());
}

BOOST_AUTO_TEST_CASE(copyMerge)
{
  ManifoldMerger mMerger;
  RealSpace R1(1);
  RealSpace R2(2);
  RealSpace R3(3);
  mMerger.addManifold(R1);
  mMerger.addManifold(R1);
  mMerger.addManifold(R2);
  mMerger.addManifold(R2);
  mMerger.addManifold(R3);
  CartesianProduct Res{&R1, &R2, &R3};
  BOOST_CHECK(mMerger.getManifold()->isSameTopology(Res));
  for (size_t i = 0; i < 3; ++i)
    BOOST_CHECK(mMerger.getManifold()->operator()(i).getInstanceId() ==
                Res(i).getInstanceId());
}

BOOST_AUTO_TEST_CASE(deepMerge)
{
  ManifoldMerger mMerger;
  RealSpace R1(1);
  RealSpace R2(2);
  RealSpace R3(3);
  RealSpace R4(4);
  RealSpace R5(5);
  RealSpace R6(6);

  CartesianProduct R1R2(R1, R2);
  CartesianProduct R3R4(R3, R4);
  CartesianProduct R5R6(R5, R6);

  mMerger.addManifold(R1R2);
  mMerger.addManifold(R3R4);
  mMerger.addManifold(R5R6);

  CartesianProduct Res{&R1, &R2, &R3, &R4, &R5, &R6};
  BOOST_CHECK(mMerger.getManifold()->isSameTopology(Res));
  for (size_t i = 0; i < 6; ++i)
    BOOST_CHECK(mMerger.getManifold()->operator()(i).getInstanceId() ==
                Res(i).getInstanceId());
}

BOOST_AUTO_TEST_CASE(deepCopyMerge)
{
  ManifoldMerger mMerger;
  RealSpace R1(1);
  RealSpace R2(2);
  RealSpace R3(3);
  RealSpace R4(4);
  RealSpace R5(5);
  RealSpace R6(6);

  CartesianProduct R1R2(R1, R2);
  CartesianProduct R2R3(R2, R3);
  CartesianProduct R3R4(R3, R4);
  CartesianProduct R4R5(R4, R5);
  CartesianProduct R5R6(R5, R6);

  mMerger.addManifold(R1R2);
  mMerger.addManifold(R1R2);
  mMerger.addManifold(R2R3);
  mMerger.addManifold(R3R4);
  mMerger.addManifold(R4R5);
  mMerger.addManifold(R5R6);

  CartesianProduct Res{&R1, &R2, &R3, &R4, &R5, &R6};
  BOOST_CHECK(mMerger.getManifold()->isSameTopology(Res));
  for (size_t i = 0; i < 6; ++i)
    BOOST_CHECK(mMerger.getManifold()->operator()(i).getInstanceId() ==
                Res(i).getInstanceId());
}
