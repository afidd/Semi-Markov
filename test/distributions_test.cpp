// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE 
//                       Agricultural Research Service 
//                  United States Department of Agriculture 
//
//   This software/database is a "United States Government Work" under the
//   terms of the United States Copyright Act.  It was written as part of
//   the author's official duties as a United States Government employee
//   and thus cannot be copyrighted.  This software/database is freely
//   available to the public for use. The Department of Agriculture (USDA)
//   and the U.S. Government have not placed any restriction on its use or
//   reproduction.
//
//   Although all reasonable efforts have been taken to ensure the accuracy
//   and reliability of the software and data, the USDA and the U.S.
//   Government do not and cannot warrant the performance or results that
//   may be obtained by using this software or data. The USDA and the U.S.
//   Government disclaim all warranties, express or implied, including
//   warranties of performance, merchantability or fitness for any
//   particular purpose.
//
//   Please cite the author in any work or product based on this material.
//
//   =========================================================================
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

  double eps=1e-9;
  GammaDistribution<RandGen> dist1(3.969, 1/1.107, 10);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(
      dist1.HazardIntegral(10,12),10)==12);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(
      dist1.HazardIntegral(11,12),11)==12);
  double hazint=dist1.HazardIntegral(11.5,12);
  double invint=dist1.ImplicitHazardIntegral(hazint, 11.5);
  BOOST_LOG_TRIVIAL(debug)<<"int "<<hazint<<" inv "<<invint;
  BOOST_CHECK(abs(12-invint)<eps);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(
      dist1.HazardIntegral(11,12.99),11)==12.99);
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



BOOST_AUTO_TEST_CASE( uniform )
{
  afidd::LogInit("debug");
  RandGen rng(1);
  UniformDistribution<RandGen> dist1(1, 3, 10);
  BOOST_CHECK(dist1.HazardIntegral(10, 10.5)==0.0);
  BOOST_CHECK(dist1.HazardIntegral(10.5, 11)==0.0);
  BOOST_CHECK(dist1.HazardIntegral(11, 11.5)>0.0);
  BOOST_CHECK(dist1.HazardIntegral(11.5, 11.99)>0.0);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(0.6, 10)>11.0);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(0.6, 11)>11.0);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(999, 11.5)>11.9);
  BOOST_LOG_TRIVIAL(debug)<<"Uniform invhazard 999, 11.5 "
      << dist1.ImplicitHazardIntegral(999, 11.5);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(999, 11.5)<=13.0);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(
      dist1.HazardIntegral(10,12),10)==12);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(
      dist1.HazardIntegral(11,12),11)==12);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(
      dist1.HazardIntegral(11.5,12),11.5)==12);
  BOOST_CHECK(dist1.ImplicitHazardIntegral(
      dist1.HazardIntegral(11,12.99),11)==12.99);
}


BOOST_AUTO_TEST_SUITE_END()
