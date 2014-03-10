#ifndef _EXPLICIT_TRANSITIONS_H_
#define _EXPLICIT_TRANSITIONS_H_ 1

#include <map>
#include <algorithm>
#include <iterator>
#include <vector>
#include <memory>
#include <type_traits>
#include "boost/mpl/for_each.hpp"
#include "boost/mpl/next_prior.hpp"
#include "logging.h"


namespace afidd
{
namespace smv
{

template<typename LM, typename State, typename RNG>
class ExplicitTransition
{
public:

  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<RNG>>>
  enabled(const State& s, const LM& lm) const
  {
    BOOST_LOG_TRIVIAL(debug) << "The base enabled is unlikely correct to call";
    return {false,
      std::unique_ptr<NoDistribution<RNG>>(new NoDistribution<RNG>())};
  }



  virtual void fire(State& s, LM& lm, RNG& rng) const
  {
    //static const auto layer_cnt=boost::mpl::size<typename LM::layers>::value;
    //static_assert(std::is_same<size_t,boost::mpl::size<typename LM::layers>::value>::value,
    //  "It's not really a size_t, is it?");
    lm.fire_by_stochiometric_coefficient(rng);
  }
};




/*! This is one way to provide a set of transitions for a GSPN.
 */
template<typename LocalMarking, typename State, typename PN, typename Random>
class ExplicitTransitions
{
public:
  typedef Random RNG;
  typedef PN PetriNet;
  typedef trans_t<PN> TransId;

  std::map<trans_t<PN>,
      std::unique_ptr<ExplicitTransition<LocalMarking,State,RNG>>>
      transitions;
};



template<typename Transitions, typename... Args>
std::pair<bool,std::unique_ptr<TransitionDistribution<typename Transitions::RNG>>>
enabled(const Transitions& et, typename Transitions::TransId trans_id,
  Args&&... args)
{
  return et.transitions.at(trans_id)->enabled(std::forward<Args>(args)...);
}



template<typename Transitions, typename... Args>
void
fire(Transitions& et, typename Transitions::TransId trans_id,
  Args&&... args)
{
  return et.transitions.at(trans_id)->fire(std::forward<Args>(args)...);
}




class EnablingPolicyNone
{
};



class EnablingPolicyInTokens
{
};


class EnablingPolicyInNoOutTokens
{
};



struct FiringPolicyInToOutHelper
{
  template<typename U>
  void operator()(U& x)
  {

  }
};


/*
class FiringPolicyInToOut
{
  void fire(const Transition<LocalMarking,State,RNG>& tr,
    State& s, LocalMarking& lm, RNG& rng)
  {
    
    for_each<typename LocalMarking::layers>(
      [&lm] (size_t layer_idx) {
        std::vector<Token> stack;
        for (auto place_idx : lm.in_places())
        {
          stack.push_back(lm.stochiometric_coefficient(place_idx));
        }
        for (auto out_idx : lm.out_places())
        {
          auto needs=-lm.stochiometric_coefficient(out_idx);
          while (needs>0)
          {
            auto available=stack.pop_back();
            auto taken=std::max(needs, available);
            lm.move<layer_idx>()
            needs-=taken;
          }
        }
      });
  }
};
*/


/*

template<typename Enable, typename Fire,
    typename LocalMarking, typename State, typename RNG>
class PolicyTransition : ExplicitTransition<LocalMarking,State,RNG>
{
public:
  virtual std::pair<bool,std::unique_ptr<Dist>>
  enabled(const Transition<LocalMarking,State,RNG>& tr, const State& s,
    const LocalMarking& lm)
  {
  }

  virtual void fire(const Transition<LocalMarking,State,RNG>& tr,
    State& s, LocalMarking& lm, RNG& rng)
  {
  }
};




template<typename Enable, typename Fire,
    typename LocalMarking, typename State, typename RNG>
class ColoredPolicyTransition : PolicyTransition<LocalMarking,State,RNG>
{
  Color _color;

public:
  ColoredPolicyTransition(const Color& c) : _color(c) {}

  virtual std::pair<bool,std::unique_ptr<Dist>>
  enabled(const Transition<LocalMarking,State,RNG>& tr, const State& s,
    const LocalMarking& lm)
  {
  }

  virtual void fire(const Transition<LocalMarking,State,RNG>& tr,
    State& s, LocalMarking& lm, RNG& rng)
  {
  }
};

*/


}
}
#endif // _EXPLICIT_TRANSITIONS_H_
