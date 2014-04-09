#ifndef _CONTINUOUS_STATE_H_
#define _CONTINUOUS_STATE_H_ 1

#include <map>
#include "stochnet.hpp"
#include "gspn.hpp"


namespace afidd
{
namespace smv
{

/*! A strategy for accumulating time intervals using just a double.
 *  This will easily underflow if there are different scales of time step.
 */
class SimpleTime
{
  double _absolute_time;

public:
  SimpleTime() {}
  SimpleTime(double start_time) : _absolute_time(start_time) {}
  double CurrentTime() const { return _absolute_time; }
  double AddTime(double interval) {
    _absolute_time+=interval;
    return _absolute_time;
  }
};




/*! Strategy to accumulate time using the Kahan summation algorithm.
 *  This is meant to ameliorate loss of precision when adding small
 *  values to larger ones.
 */
class KahanTime
{
  double _sum;
  double _compensation;

public:
  KahanTime() : _sum(0), _compensation(0) {}
  KahanTime(double init) : _sum(init), _compensation(0) {}
  double CurrentTime() const { return _sum; }
  double AddTime(double interval) {
    auto adjusted_interval=interval - _compensation;
    auto next_sum=_sum + adjusted_interval;
    _compensation=(next_sum - _sum) - adjusted_interval;
    _sum=next_sum;
    return _sum;
  }
};




template<typename Mark, typename UserState=detail::NoExtraState,
    typename TimeStrategy=KahanTime>
class GSPNState : public TimeStrategy
{
public:
  typedef Mark Marking;

  Marking marking;
  UserState user;
};

}
}



#endif // _CONTINUOUS_STATE_H_
