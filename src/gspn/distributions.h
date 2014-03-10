#ifndef _DISTRIBUTIONS_H_
#define _DISTRIBUTIONS_H_ 1

#include <tuple>
#include <limits>
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

} // smv
} // afidd
#endif // _DISTRIBUTIONS_H_
