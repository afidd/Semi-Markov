#include <type_traits>
#include <utility>
#include <algorithm>
#include <random>
#include <sstream>
#include <fstream>
#define BOOST_TEST_DYN_LINK 1
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "logging.hpp"
#include "distributions.hpp"


using namespace afidd::smv;

using RandGen=std::mt19937;


BOOST_AUTO_TEST_SUITE( distributions )

BOOST_AUTO_TEST_CASE( exponential )
{
  afidd::LogInit("debug");
  RandGen rng(1);
  for (auto lambda : std::vector<double>{1.0, 20.0, 0.01}) {
    ExponentialDistribution<RandGen> dist1(lambda, 0.0);
    std::vector<double> vals(100000);
    for (auto it=vals.begin(); it!=vals.end(); ++it) {
      *it=dist1.Sample(0.0, rng);
    }
    BOOST_CHECK(dist1.CheckSamples(vals, 0.0));
  }

}



BOOST_AUTO_TEST_CASE( weibull )
{
  afidd::LogInit("debug");
  RandGen rng(1);
  for (auto lambda : std::vector<double>{1.0, 20.0, 0.01}) {
    WeibullDistribution<RandGen> dist1(lambda, 1.0, 0.0);
    std::vector<double> vals(100000);
    for (auto it=vals.begin(); it!=vals.end(); ++it) {
      *it=dist1.Sample(0.0, rng);
    }

    std::sort(vals.begin(), vals.end());

    std::stringstream name;
    name<<"weibull_"<<lambda<<".txt";
    std::ofstream output(name.str());
    size_t cnt=0;
    output<<0.0<<" "<<0.0<<std::endl;
    for (auto v : vals) {
      ++cnt;
      output<<v<<" "<<cnt/((double)vals.size())<<std::endl;
    }
    BOOST_CHECK(dist1.CheckSamples(vals, 0.0));
  }
}



BOOST_AUTO_TEST_CASE( gamma )
{
  afidd::LogInit("debug");
  RandGen rng(1);
  for (auto lambda : std::vector<double>{1.0, 2.0, 10.0}) {
    GammaDistribution<RandGen> dist1(lambda, 0.5, 0.0);
    std::vector<double> vals(100000);
    for (auto it=vals.begin(); it!=vals.end(); ++it) {
      *it=dist1.Sample(0.0, rng);
    }

    std::sort(vals.begin(), vals.end());

    std::stringstream name;
    name<<"gamma_"<<lambda<<".txt";
    std::ofstream output(name.str());
    size_t cnt=0;
    output<<0.0<<" "<<0.0<<std::endl;
    for (auto v : vals) {
      ++cnt;
      output<<v<<" "<<cnt/((double)vals.size())<<std::endl;
    }
    BOOST_CHECK(dist1.CheckSamples(vals, 0.0));
  }
}


void tofile(const std::vector<double>& vals, const std::string& name)
{
  std::ofstream output(name);
  size_t cnt=0;
  output<<0.0<<" "<<0.0<<std::endl;
  for (auto v : vals) {
    ++cnt;
    output<<v<<" "<<cnt/((double)vals.size())<<std::endl;
  }
}




BOOST_AUTO_TEST_CASE( piecewiselinear )
{
  afidd::LogInit("debug");
  RandGen rng(1);
  std::vector<double> b={0,   1,   2, 3, 3.2, 5};
  std::vector<double> w={0.2, 1, 0.5, 0,   0, 1};

  double cutoff=2;
  PiecewiseLinearDistribution<RandGen> dist1(b, w, 0.0);
  std::vector<double> vals(10000);
  for (auto it=vals.begin(); it!=vals.end(); ++it) {
    *it=dist1.Sample(cutoff, rng);
  }

  std::sort(vals.begin(), vals.end());
  tofile(vals, "piecewise.txt");
  BOOST_CHECK(dist1.CheckSamples(vals, 0.0));

  std::piecewise_linear_distribution<double> pld(b.begin(), b.end(), w.begin());
  for (auto pt=vals.begin(); pt!=vals.end(); ++pt) {
    double possible;
    for (possible=pld(rng); possible<cutoff; possible=pld(rng))
    {}
    *pt=possible;
  }
  std::sort(vals.begin(), vals.end());
  tofile(vals, "stdpiece.txt");
}


BOOST_AUTO_TEST_SUITE_END()
