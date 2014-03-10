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
template<typename LocalMarking, typename State,
    typename PetriGraph, typename Random>
class ExplicitTransitions
{
public:
  typedef Random RNG;
  typedef trans_t<PetriGraph> TransId;
  ExplicitTransitions(PetriGraph& graph) : graph(std::move(graph)) {}

  std::map<trans_t<PetriGraph>,
      std::unique_ptr<ExplicitTransition<LocalMarking,State,RNG>>>
      transitions;
  PetriGraph graph;
};


template<typename LocalMarking, typename State,
    typename PetriGraph, typename Random>
struct petri_place<ExplicitTransitions<LocalMarking,State,PetriGraph,Random>>
{
  typedef typename boost::graph_traits<PetriGraph>::vertex_descriptor type;
};



template<typename LocalMarking, typename State,
    typename PetriGraph, typename Random>
struct petri_transition<ExplicitTransitions<LocalMarking,State,PetriGraph,Random>>
{
  typedef typename boost::graph_traits<PetriGraph>::vertex_descriptor type;
};



template<typename LocalMarking, typename State,
    typename PetriGraph, typename Random>
std::vector<std::tuple<size_t,size_t,int>>
neighbors_of_transition(
  ExplicitTransitions<LocalMarking,State,PetriGraph,Random>& et,
  typename petri_transition<PetriGraph>::type trans_id)
{
  return neighbors_of_transition(et.graph, trans_id);
}


template<typename F, typename LocalMarking, typename State,
    typename PetriGraph, typename Random>
void neighbors_of_places(
  ExplicitTransitions<LocalMarking,State,PetriGraph,Random>& et,
  const std::set<typename petri_place<PetriGraph>::type>& place_id, F func)
{
  return neighbors_of_places(et.graph, place_id, func);
}


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


}
}
#endif // _EXPLICIT_TRANSITIONS_H_
