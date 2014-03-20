#include <type_traits>
#include <utility>
#include <random>
#define BOOST_TEST_DYN_LINK 1
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "logging.h"
#include "distributions.h"


using namespace afidd::smv;

using RandGen=std::mt19937;


BOOST_AUTO_TEST_SUITE( distributions )

BOOST_AUTO_TEST_CASE( exponential )
{
  afidd::log_init("debug");
  RandGen rng(1);
  for (auto lambda : std::vector<double>{1.0, 20.0, 0.01})
  {
    ExponentialDistribution<RandGen> dist1(lambda);
    std::vector<double> vals(1000000);
    for (auto it=vals.begin(); it!=vals.end(); ++it)
    {
      *it=dist1.sample(0.0, 0.0, rng);
    }
    BOOST_CHECK(dist1.check_samples(vals, 0.0));
  }

}



BOOST_AUTO_TEST_CASE( weibull )
{
  afidd::log_init("debug");
  RandGen rng(1);
  for (auto lambda : std::vector<double>{1.0, 20.0, 0.01})
  {
    WeibullDistribution<RandGen> dist1(lambda, 1.0);
    std::vector<double> vals(1000000);
    for (auto it=vals.begin(); it!=vals.end(); ++it)
    {
      *it=dist1.sample(0.0, 0.0, rng);
    }
    BOOST_CHECK(dist1.check_samples(vals, 0.0));
  }

}



BOOST_AUTO_TEST_SUITE_END()
