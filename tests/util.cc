// Copyright (C) 2015 by Benjamin Chrétien, CNRS-LIRMM.
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

#include <memory>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <roboptim/core/manifold-map/util.hh>

boost::shared_ptr<boost::test_tools::output_test_stream> output;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (util)
{
  using namespace roboptim;

  output = retrievePattern ("util");

  boost::shared_ptr<int> a = boost::make_shared<int> (42);
  BOOST_CHECK_EQUAL (a.use_count (), 1);

  // Note: ref counters are not shared.
  std::shared_ptr<int> b = toStd (a);
  BOOST_CHECK_EQUAL (a.use_count (), 2);
  BOOST_CHECK_EQUAL (b.use_count (), 1);

  boost::shared_ptr<int> c = toBoost (b);
  BOOST_CHECK_EQUAL (a.use_count (), 2);
  BOOST_CHECK_EQUAL (b.use_count (), 2);
  BOOST_CHECK_EQUAL (c.use_count (), 1);

  (*output) << *a << std::endl;
  (*output) << *b << std::endl;
  (*output) << *c << std::endl;

  // Change c
  *c = 12;
  (*output) << *a << std::endl;
  (*output) << *b << std::endl;
  (*output) << *c << std::endl;

  // Delete b
  b.reset ();
  BOOST_CHECK_EQUAL (a.use_count (), 2);
  BOOST_CHECK_EQUAL (c.use_count (), 1);
  (*output) << *a << std::endl;
  (*output) << *c << std::endl;

  int foo = 1337;
  (*output) << foo << std::endl;

  // If you actually do this, you will burn in a very special level of
  // hell...
  std::shared_ptr<int> bar (&foo, NoopDeleter<int> ());
  (*output) << *bar << std::endl;
  bar.reset ();
  (*output) << foo << std::endl;

  std::cout << output->str () << std::endl;
  BOOST_CHECK (output->match_pattern ());
}

BOOST_AUTO_TEST_SUITE_END ()
