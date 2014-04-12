#ifndef _CONTINUOUS_STATE_H_
#define _CONTINUOUS_STATE_H_ 1

#include <map>
#include "stochnet.hpp"
#include "gspn.hpp"


namespace afidd
{
namespace smv
{
template<typename Mark, typename TransitionKey,
    typename UserState=detail::NoExtraState>
class GSPNState
{
public:
  typedef Mark Marking;
  Marking marking;
  double absolute_time;
  double CurrentTime() const { return absolute_time; }
  double SetTime(double timeprime) { absolute_time=timeprime; }
  UserState user;
  TransitionKey last_transition;
};

}
}



#endif // _CONTINUOUS_STATE_H_
