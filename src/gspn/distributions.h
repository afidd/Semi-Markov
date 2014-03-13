#ifndef _DISTRIBUTIONS_H_
#define _DISTRIBUTIONS_H_ 1

#include <tuple>
#include <limits>
#include "boost/math/distributions/gamma.hpp"
#include "stochnet.h"
#include "gspn_random.h"


namespace afidd
{
namespace smv
{



template<typename RNG>
class TransitionDistribution
{
public:
  virtual double sample(double enabling_time, double current_time,
      RNG& rng) const=0;
};



template<typename RNG>
class NoDistribution : public TransitionDistribution<RNG>
{
public:
  virtual double sample(double enabling_time, double current_time,
      RNG& rng) const { return std::numeric_limits<double>::infinity(); };
};


template<typename RNG>
class ExponentialDistribution : public TransitionDistribution<RNG>
{
  using ParamType=std::tuple<double>;
  ParamType _params;

public:
  ExponentialDistribution(double lambda)
  : _params(lambda) {}

  virtual double sample(double enabling_time,
      double current_time, RNG& rng) const
  {
    return -std::get<0>(_params)*std::log(uniform(rng));
  }

  double sample_vector(
    const std::vector<ParamType> params,
    const std::vector<double>& enabling_time,
    double current_time, RNG& rng)
  {

  }

};



template<typename RNG>
class WeibullDistribution : public TransitionDistribution<RNG>
{
  using ParamType=std::tuple<double,double>;
  ParamType _params;
public:
  WeibullDistribution(double lambda, double k) : _params{lambda, k} {}

  virtual double sample(double enabling_time, double current_time,
      RNG& rng) const
  {
    double l=std::get<0>(_params);
    double k=std::get<1>(_params);
    double d=current_time-enabling_time;
    double U=uniform(rng);

    if (d>0)
    {
      return l*std::pow(-std::log(1-U)+std::pow(d/l, k), 1/k)-d;
    }
    else
    {
      return l*std::pow(-std::log(1-U), 1/k);
    }
  };
};




template<typename RNG>
class GammaDistribution : public TransitionDistribution<RNG>
{
  using Params=std::tuple<double,double>;
  Params _params;
public:
  GammaDistribution(double alpha, double beta)
  : _params{alpha, beta}
  {
  }


  virtual double sample(double enabling_time, double current_time,
      RNG& rng) const
  {
    double d=current_time-enabling_time;
    double U=uniform(rng);
    auto dist=boost::math::gamma_distribution<double>(
      std::get<0>(_params), std::get<1>(_params));

    if (d>0)
    {
      return boost::math::quantile(dist, U);
    }
    else
    {
      auto cumulative=boost::math::cdf(dist, d);
      return boost::math::quantile(dist, U*(1-cumulative) + cumulative) - d;
    }
  };
};


} // smv
} // afidd
#endif // _DISTRIBUTIONS_H_
