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
#ifndef _DISTRIBUTIONS_H_
#define _DISTRIBUTIONS_H_ 1

#include <tuple>
#include <limits>
#include "boost/math/distributions/gamma.hpp"
#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/mean.hpp"
#include "boost/accumulators/statistics/moment.hpp"
#include "boost/accumulators/statistics/variance.hpp"
#include "logging.hpp"
#include "stochnet.hpp"
#include "gspn_random.hpp"


namespace afidd
{
namespace smv
{


namespace detail
{
  double FracError(double a, double b) {
    return std::abs((a-b)/a);
  }

  bool CheckFracError(double a, double b, double tol, const std::string& m) {
    if (detail::FracError(a, b)>tol)
    {
      BOOST_LOG_TRIVIAL(info) << "Fractional error of "<<m<< " too large. "
        "Expected "<< a << " found "<<b;
        return false;
    }
    return true;
  }
}


template<typename RNG>
class TransitionDistribution
{
public:
  virtual double Sample(double current_time, RNG& rng) const=0;
  virtual double EnablingTime() const=0;
};



template<typename RNG>
class NoDistribution : public TransitionDistribution<RNG>
{
public:
  virtual double Sample(double current_time, RNG& rng) const {
    return std::numeric_limits<double>::infinity(); };
  virtual double EnablingTime() const { return 0; }
};



template<typename RNG>
class ExponentialDistribution : public TransitionDistribution<RNG>
{
  using ParamType=std::tuple<double, double, double>;
  ParamType params_;

public:
  ExponentialDistribution(double lambda, double enabling_time,
    double normal=1.0) : params_(lambda, enabling_time, normal)
  {}



  virtual double Sample(double current_time, RNG& rng) const {
    auto U=uniform(rng)/std::get<2>(params_);
    if (U>=1) {
      return std::numeric_limits<double>::infinity();
    }
    return -std::log(U)/std::get<0>(params_);
  }


  virtual double EnablingTime() const {
    return std::get<1>(params_);
  }



  /*! Check whether generated samples fit expectations.
   *  dt is the difference between enabling_time and current_time,
   *  assumed constant for all samples.
   */
  bool CheckSamples(const std::vector<double>& samples, double dt) {
    namespace bac=boost::accumulators;
    bool pass=true;
    bac::accumulator_set<double, bac::stats<bac::tag::mean,
        bac::tag::moment<2>, bac::tag::moment<3>,
        bac::tag::variance(bac::lazy)>> acc;
    for (auto s : samples) {
      acc(s);
    }

    double lambda_estimator=1/bac::mean(acc);
    double lambda=std::get<0>(params_);
    auto too_low=lambda < lambda_estimator*(1-1.96/std::sqrt(samples.size()));
    auto too_high=lambda > lambda_estimator*(1+1.96/std::sqrt(samples.size()));
    if (too_low || too_high) {
      BOOST_LOG_TRIVIAL(info)<<"Parameter not in bounds. Low? "<<
        too_low << " high? " << too_high;
        pass=false;
    }

    double variance=bac::variance(acc);
    // An estimator for the variance?
    pass=detail::CheckFracError(variance, std::pow(lambda, -2), 0.01,
      "variance");
    return pass;
  }
};






template<typename RNG>
class ShiftedExponentialDistribution : public TransitionDistribution<RNG>
{
  using ParamType=std::tuple<double, double, double, double>;
  ParamType params_;

public:
  ShiftedExponentialDistribution(double lambda, double enabling_time,
    double shift=0.0, double normal=1.0)
  : params_(lambda, enabling_time, shift, normal) {}



  virtual double Sample(double current_time, RNG& rng) const {
    auto U=uniform(rng)/std::get<3>(params_);
    double te=std::get<1>(params_);
    double ts=std::get<2>(params_);
    if (U>=1) {
      return std::numeric_limits<double>::infinity();
    }
    if (current_time>te+ts) {
      return -std::log(U)/std::get<0>(params_);
    } else {
      return te+ts-std::log(U)/std::get<0>(params_)-current_time;
    }
  }


  virtual double EnablingTime() const {
    return std::get<1>(params_);
  }


  /*! Check whether generated samples fit expectations.
   *  dt is the difference between enabling_time and current_time,
   *  assumed constant for all samples.
   */
  bool CheckSamples(const std::vector<double>& samples, double dt) {
    namespace bac=boost::accumulators;
    bool pass=true;
    bac::accumulator_set<double, bac::stats<bac::tag::mean,
        bac::tag::moment<2>, bac::tag::moment<3>,
        bac::tag::variance(bac::lazy)>> acc;
    for (auto s : samples) {
      acc(s);
    }

    double lambda_estimator=1/bac::mean(acc);
    double lambda=std::get<0>(params_);
    auto too_low=lambda < lambda_estimator*(1-1.96/std::sqrt(samples.size()));
    auto too_high=lambda > lambda_estimator*(1+1.96/std::sqrt(samples.size()));
    if (too_low || too_high) {
      BOOST_LOG_TRIVIAL(info)<<"Parameter not in bounds. Low? "<<
        too_low << " high? " << too_high;
        pass=false;
    }

    double variance=bac::variance(acc);
    // An estimator for the variance?
    pass=detail::CheckFracError(variance, std::pow(lambda, -2), 0.01,
      "variance");
    return pass;
  }

};






template<typename RNG>
class WeibullDistribution : public TransitionDistribution<RNG>
{
  using ParamType=std::tuple<double,double,double,double,double>;
  ParamType params_;
public:
  WeibullDistribution(double lambda, double k, double enabling_time,
    double shift=0.0, double normal=1.0)
  : params_{lambda, k, enabling_time, shift, normal} {}

  virtual double Sample(double current_time,
      RNG& rng) const
  {
    double l=std::get<0>(params_);
    double k=std::get<1>(params_);
    double enabling_time=std::get<2>(params_);
    double shift=std::get<3>(params_);
    double U=uniform(rng)/std::get<4>(params_);
    if (U>=1) {
      return std::numeric_limits<double>::infinity();
    }

    double d=current_time-(enabling_time+shift);
    if (d>0) {
      return l*std::pow(-std::log(1-U)+std::pow(d/l, k), 1/k)-d;
    } else {
      return -d+l*std::pow(-std::log(1-U), 1/k);
    }
  }



  virtual double EnablingTime() const {
    return std::get<2>(params_);
  }


  /*! Check estimators to see if the distribution is correct.
   *  \lambda^k = \frac{1}{N}\sum_{i=1}^N(x_i^k-x_N^k)
   *  where x_N is the largest observed sample.
   *  From wikipedia. Yup.
   */
  bool CheckSamples(const std::vector<double>& samples, double dt) {
    namespace bac=boost::accumulators;
    bool pass=true;
    double l=std::get<0>(params_);
    double k=std::get<1>(params_);


    bac::accumulator_set<double, bac::stats<bac::tag::mean,
        bac::tag::moment<2>, bac::tag::moment<3>,
        bac::tag::variance(bac::lazy)>> acc;
    for (auto st : samples) {
      acc(st);
    }

    if (std::abs(dt-0)<0.0000001) {
      double expected_mean=l*std::tgamma(1+1/k);
      double mean=bac::mean(acc);
      if (detail::FracError(expected_mean, mean)>0.01) {
        pass=false;
        BOOST_LOG_TRIVIAL(info)<<"Expected mean " << expected_mean
          << " but found "<< mean;
      }

      double expected_variance=l*l*(
        std::tgamma(1+2/k)-std::pow(std::tgamma(1+1/k), 2)
        );
      double variance=bac::variance(acc);
      if (detail::FracError(expected_variance, variance)>0.01) {
        pass=false;
        BOOST_LOG_TRIVIAL(info)<<"Expected variance " << expected_variance
          << " but found "<< variance;
      }

      double min=*std::min_element(samples.begin(), samples.end());
      double mink=std::pow(min, k);

      double total=0.0;
      for (auto s : samples) {
        total+=std::pow(s, k);
      }
      double l_estimator=std::pow(total/samples.size() - mink, 1/k);

      if (detail::FracError(l, l_estimator) > 0.01) {
        pass=false;
        BOOST_LOG_TRIVIAL(info)<<"Expected lambda " << l
          << " but found "<<l_estimator;
      }

      double numerator=0.0;
      double denominator=0.0;
      double logsum=0.0;
      for (auto t : samples) {
        numerator+=std::pow(t, k)*std::log(t) - mink*std::log(min);
        denominator+=std::pow(t, k)-mink;
        logsum+=std::log(t);
      }
      double k_est_inv=numerator/denominator - logsum/samples.size();
      double k_est=1.0/k_est_inv;

      pass=detail::CheckFracError(k, k_est, 0.01, "k");
    }

    return pass;
  }
};




template<typename RNG>
class GammaDistribution : public TransitionDistribution<RNG>
{
  using Params=std::tuple<double,double,double, double, double>;
  Params params_;
public:
  GammaDistribution(double alpha, double theta, double enabling_time,
      double shift=0.0, double normal=1.0)
  : params_{alpha, theta, enabling_time, shift, normal}
  {}


  virtual double Sample(double current_time, RNG& rng) const {
    double te=std::get<2>(params_);
    double ts=std::get<3>(params_);
    double t0=current_time;
    double U=uniform(rng)/std::get<4>(params_);
    auto dist=boost::math::gamma_distribution<double>(
      std::get<0>(params_), std::get<1>(params_));
    if (U>=1) {
      return std::numeric_limits<double>::infinity();
    }

    double d=t0-(te+ts);
    if (d>0) {
      auto cumulative=boost::math::cdf(dist, d);
      return boost::math::quantile(dist, U*(1-cumulative) + cumulative) - d;
    } else {
      return boost::math::quantile(dist, U)-d;
    }
  }



  virtual double EnablingTime() const {
    return std::get<2>(params_);
  }



  bool CheckSamples(const std::vector<double>& samples, double dt) {
    namespace bac=boost::accumulators;
    bool pass=true;
    double a=std::get<0>(params_);
    double th=std::get<1>(params_);

    bac::accumulator_set<double, bac::stats<bac::tag::mean,
        bac::tag::moment<2>, bac::tag::moment<3>,
        bac::tag::variance(bac::lazy)>> acc;
    for (auto st : samples) {
      acc(st);
    }

    if (std::abs(dt-0)<0.0000001) {
      double expected_mean=a*th;
      double expected_variance=a*th*th;
      double expected_skew=2/std::sqrt(a);
      pass=detail::CheckFracError(
          expected_mean, bac::mean(acc), 0.01, "mean");
      pass=detail::CheckFracError(
          expected_variance, bac::variance(acc), 0.01, "variance");
      pass=detail::CheckFracError(
          expected_skew, bac::moment<3>(acc), 0.01, "skew");
    }

    double th_est=bac::mean(acc)/a;
    pass=detail::CheckFracError(th, th_est, 0.01, "theta");

    // Following wikipedia on Gamma distribution...
    double slog=0.0;
    for (auto sl : samples) {
      slog+=std::log(sl);
    }
    double s=std::log(bac::mean(acc))-slog/samples.size();
    double a_est=(3-s+std::sqrt((s-3)*(s-3)+24*s))/(12*s);
    pass=detail::CheckFracError(a, a_est, 0.03, "alpha");
    return pass;
  }
};





/*! A sequence of contiguous subintervals defined by endpoints, b,
 *  and height at each point, w. The distribution is zero outside
 *  this range.
 */
template<typename RNG>
class PiecewiseLinearDistribution : public TransitionDistribution<RNG>
{
  using Params=std::tuple<std::vector<double>,std::vector<double>,double,
      double,double>;
  Params params_;

  std::vector<double> _partial_sum;

public:
  PiecewiseLinearDistribution(const std::vector<double>& b,
    const std::vector<double>& w, double enabling_time,
    double shift=0.0, double normal=1.0)
  : params_{b, w, enabling_time, shift, normal}, _partial_sum(b.size()) {
    assert(b.size()>1);
    assert(b.size()==w.size());
    assert(std::is_sorted(b.begin(), b.end()));

    double total=0;
    for (int idx=0; idx<b.size()-1; ++idx) {
      _partial_sum[idx]=total;
      total+=0.5*(w[idx]+w[idx+1])*(b[idx+1]-b[idx]);
    }
    _partial_sum[b.size()-1]=total;
  }



  virtual double EnablingTime() const {
    return std::get<2>(params_);
  }


  virtual double Sample(double current_time, RNG& rng) const {
    const auto& b=std::get<0>(params_);
    const auto& w=std::get<1>(params_);
    double enabling_time=std::get<2>(params_);
    double from_time=current_time-enabling_time;
    double U=uniform(rng)/std::get<4>(params_);

    // What is the first b not less than current time? Store as b_{i+1}.
    double S=1.0;
    double within_integral=0.0;
    auto ip1_iter=std::lower_bound(b.begin(), b.end(), from_time);
    typename std::iterator_traits<decltype(ip1_iter)>::difference_type p0{0};
    if (ip1_iter==b.begin()) {
      // time is before start of piecewise. Return the whole.
      S=_partial_sum[b.size()-1];
    } else if (ip1_iter==b.end()) {
      BOOST_LOG_TRIVIAL(warning) << "Normalization of piecewise distribution "
          "indicates that it should have fired.";
      return 0;
    } else {
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
    if (bp1_iter==_partial_sum.end()) {
      BOOST_LOG_TRIVIAL(info)<<"Piecewise distribution over bounds.";
      return std::numeric_limits<double>::infinity();
    } else {
      // Solve the quadratic for the second root and then
      // substitute back into r+ x r- = c/a. This way we can
      // handle zero slope and zero weight.
      auto j=(bp1_iter-_partial_sum.begin())-1;
      c-=_partial_sum[j];
      double m=(w[j+1]-w[j])/(b[j+1]-b[j]); // can be zero.
      double t=b[j]+2*c/(w[j]+std::sqrt(w[j]*w[j]+2*m*c));
      if (t<=b[b.size()-1] && t>=0) {
        return t;
      } else {
        BOOST_LOG_TRIVIAL(error)<<"infinite piecewise j="<<j
          <<" m="<<m<<" c="<<c<<" w="<<w[j]<<" w+1="<<w[j+1]
          <<" discr="<<w[j]*w[j]+2*m*c<<" t="<<t;
        return 0;
      }
    }
  }


  double CumulativeProbability(double enabling_time, double current_time,
    double x)
  {
    return 0.0;
  }


  bool CheckSamples(const std::vector<double>& samples, double dt)
  {
    bool pass=true;

    // Try Kolmogorov-Smirnov test.

    return pass;
  }
};


} // smv
} // afidd
#endif // _DISTRIBUTIONS_H_
