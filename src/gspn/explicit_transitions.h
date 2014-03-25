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
  enabled(const State& s, const LM& lm, double current_time) const
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
template<typename ETState,typename Random>
class ExplicitTransitions
{
  // The GSPN gets the marking type from the State.
  typedef typename ETState::Marking ETMarking;
  typedef PetriGraphType PetriGraph;

public:
  typedef smv::LocalMarking<ETMarking> LocalMarking;
  typedef Random RNG;
  typedef boost::graph_traits<PetriGraph>::vertex_descriptor place_type;
  typedef boost::graph_traits<PetriGraph>::vertex_descriptor transition_type;
  typedef smv::ExplicitTransition<LocalMarking,ETState,RNG> Transition;

  ExplicitTransitions(PetriGraph& graph) : graph(std::move(graph)) {}
  ExplicitTransitions(const ExplicitTransitions&)=delete;
  ExplicitTransitions& operator=(const ExplicitTransitions&)=delete;


  ExplicitTransitions(ExplicitTransitions&& other)
  {
    transitions=std::move(other.transitions);
    graph=std::move(other.graph);
  }


  ExplicitTransitions& operator=(ExplicitTransitions&& other)
  {
    if (this!=other)
    {
      transitions=std::move(other.transitions);
      graph=std::move(other.graph);
    }
  }

  ~ExplicitTransitions() {}

  std::map<transition_type,
      std::unique_ptr<ExplicitTransition<LocalMarking,ETState,RNG>>>
      transitions;
  PetriGraph graph;
};


template<typename State, typename Random>
struct petri_place<ExplicitTransitions<State,Random>>
{
  typedef typename ExplicitTransitions<State,Random>::place_type type;
};



template<typename State, typename Random>
struct petri_transition<ExplicitTransitions<State,Random>>
{
  typedef typename ExplicitTransitions<State,Random>::transition_type type;
};



template<typename State, typename Random>
std::vector<std::tuple<size_t,size_t,int>>
neighbors_of_transition(
  ExplicitTransitions<State,Random>& et,
  typename ExplicitTransitions<State,Random>::transition_type trans_id)
{
  return neighbors_of_transition(et.graph, trans_id);
}


template<typename F, typename State, typename Random>
void neighbors_of_places(
  ExplicitTransitions<State,Random>& et,
  const std::set<typename ExplicitTransitions<State,Random>::place_type>& place_id, F func)
{
  return neighbors_of_places(et.graph, place_id, func);
}


template<typename Transitions, typename... Args>
std::pair<bool,std::unique_ptr<TransitionDistribution<typename Transitions::RNG>>>
enabled(const Transitions& et, typename Transitions::transition_type trans_id,
  Args&&... args)
{
  return et.transitions.at(trans_id)->enabled(std::forward<Args>(args)...);
}



template<typename Transitions, typename... Args>
void
fire(Transitions& et, typename Transitions::transition_type trans_id,
  Args&&... args)
{
  et.transitions.at(trans_id)->fire(std::forward<Args>(args)...);
}


}
}
#endif // _EXPLICIT_TRANSITIONS_H_
