#ifndef _CONTINUOUS_STATE_H_
#define _CONTINUOUS_STATE_H_ 1

#include <map>
#include "stochnet.h"
#include "gspn.h"


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
  double current_time() const { return _absolute_time; }
  double add_time(double interval)
  {
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
  double current_time() const { return _sum; }
  double add_time(double interval)
  {
    auto adjusted_interval=interval - _compensation;
    auto next_sum=_sum + adjusted_interval;
    _compensation=(next_sum - _sum) - adjusted_interval;
    _sum=next_sum;
    return _sum;
  }
};




template<typename PN, typename Mark, typename TimeStrategy=KahanTime>
class State : public TimeStrategy
{
public:
  typedef Mark Marking;

  Marking marking;
  std::map<trans_t<PN>,double> enabling_time;
};

}
}



#endif // _CONTINUOUS_STATE_H_
