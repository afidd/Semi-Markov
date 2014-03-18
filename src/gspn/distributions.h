#ifndef _DISTRIBUTIONS_H_
#define _DISTRIBUTIONS_H_ 1

#include <tuple>
#include <limits>
#include "boost/math/distributions/gamma.hpp"
#include "logging.h"
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
  using ParamType=std::tuple<double, double>;
  ParamType _params;

public:
  ExponentialDistribution(double lambda, double normal=1.0)
  : _params(lambda, normal) {}

  virtual double sample(double enabling_time,
      double current_time, RNG& rng) const
  {
    auto U=uniform(rng)/std::get<1>(_params);
    if (U>=1)
    {
      return std::numeric_limits<double>::infinity();
    }
    return -std::get<0>(_params)*std::log(U);
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
  using ParamType=std::tuple<double,double,double>;
  ParamType _params;
public:
  WeibullDistribution(double lambda, double k, double normal=1.0)
  : _params{lambda, k, normal} {}

  virtual double sample(double enabling_time, double current_time,
      RNG& rng) const
  {
    double l=std::get<0>(_params);
    double k=std::get<1>(_params);
    double d=current_time-enabling_time;
    double U=uniform(rng)/std::get<2>(_params);
    if (U>=1)
    {
      return std::numeric_limits<double>::infinity();
    }

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
  using Params=std::tuple<double,double,double>;
  Params _params;
public:
  GammaDistribution(double alpha, double beta, double normal=1.0)
  : _params{alpha, beta, normal}
  {
  }


  virtual double sample(double enabling_time, double current_time,
      RNG& rng) const
  {
    double d=current_time-enabling_time;
    double U=uniform(rng)/std::get<2>(_params);
    auto dist=boost::math::gamma_distribution<double>(
      std::get<0>(_params), std::get<1>(_params));
    if (U>=1)
    {
      return std::numeric_limits<double>::infinity();
    }

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





/*! A sequence of contiguous subintervals defined by endpoints, b,
 *  and height at each point, w. The distribution is zero outside
 *  this range.
 */
template<typename RNG>
class PiecewiseLinearDistribution : public TransitionDistribution<RNG>
{
  using Params=std::tuple<std::vector<double>,std::vector<double>,double>;
  Params _params;

  std::vector<double> _partial_sum;

public:
  PiecewiseLinearDistribution(const std::vector<double>& b,
    const std::vector<double>& w, double normal=1.0)
  : _params{b, w, normal}, _partial_sum(b.size())
  {
    assert(b.size()>1);
    assert(b.size()==w.size());
    assert(std::is_sorted(b.begin(), b.end()));

    double total=0;
    for (int idx=0; idx<b.size()-2; ++idx)
    {
      _partial_sum[idx]=total;
      total+=0.5*(w[idx]+w[idx+1])*(b[idx+1]-b[idx]);
    }
    _partial_sum[b.size()-1]=total;
  }



  virtual double sample(double enabling_time, double current_time,
    RNG& rng) const
  {
    const auto& b=std::get<0>(_params);
    const auto& w=std::get<1>(_params);

    double from_time=current_time-enabling_time;
    double U=uniform(rng)/std::get<2>(_params);

    // What is the first b not less than current time? Store as b_{i+1}.
    double S=1.0;
    double within_integral=0.0;
    auto ip1_iter=std::lower_bound(b.begin(), b.end(), from_time);
    typename std::iterator_traits<decltype(ip1_iter)>::difference_type p0;
    if (ip1_iter==b.begin())
    {
      // time is before start of piecewise. Return the whole.
      S=_partial_sum[b.size()-1];
    }
    else if (ip1_iter==b.end())
    {
      BOOST_LOG_TRIVIAL(warning) << "Normalization of piecewise distribution "
          "indicates that it should have fired.";
      return 0;
    }
    else
    {
      auto p1=(ip1_iter-b.begin());
      p0=p1-1;
      double dt=from_time-b[p0];
      within_integral=_partial_sum[p0]+
          dt* (w[p0]+ 0.5*dt*(w[p1]-w[p0])/(b[p1]-b[p0]));
      S=_partial_sum[b.size()-1]-within_integral;
    }

    double c=S*U+within_integral;
    auto bp1_iter=std::lower_bound(_partial_sum.begin()+p0, _partial_sum.end(),
        c);
    if (bp1_iter==_partial_sum.end())
    {
      return std::numeric_limits<double>::infinity();
    }
    else
    {
      // Solve the quadratic for the second root and then
      // substitute back into r+ x r- = c/a. This way we can
      // handle zero slope and zero weight.
      auto j=(bp1_iter-_partial_sum.begin())-1;
      double m=(w[j+1]-w[j])/(b[j+1]-b[j]); // can be zero.
      double t=b[j]+2*c/(w[j]+std::sqrt(w[j]*w[j]+2*m*c));
      return t;
    }
  }
};


} // smv
} // afidd
#endif // _DISTRIBUTIONS_H_
