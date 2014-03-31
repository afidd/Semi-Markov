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
#include "build_graph.h"
#include "gspn.h"
#include "logging.h"


namespace afidd
{
namespace smv
{


template<typename LM, typename RNG, typename ExtraState=detail::NoExtraState>
class ExplicitTransition
{
public:
  typedef ExtraState UserState;

  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<RNG>>>
  enabled(const ExtraState& s, const LM& lm,
          double enabling_time, double current_time) const
  {
    BOOST_LOG_TRIVIAL(debug) << "The base enabled is unlikely correct to call";
    return {false,std::unique_ptr<TransitionDistribution<RNG>>(nullptr)};
  }



  virtual void fire(ExtraState& s, LM& lm, RNG& rng) const
  {
    lm.template transfer_by_stochiometric_coefficient<0>(rng);
  }
};




/*! This is one way to provide a set of transitions for a GSPN.
 *  ETState is the state of the system. PKey is a place key, the identifier
 *  for the place. TKey is the identifier for the transition. Random
 *  is the random number generator.
 */
template<typename PKey, typename TKey, typename Local,
  typename ETRand, typename ExtraState=detail::NoExtraState>
class ExplicitTransitions
{
  typedef PetriGraphType PetriGraph;
  typedef BiGraphCorrespondence<PKey,TKey,
    boost::graph_traits<PetriGraph>::vertex_descriptor> BiMap;

public:
  typedef PKey PlaceKey;
  typedef TKey TransitionKey;
  typedef ETRand RNG;
  typedef boost::graph_traits<PetriGraph>::vertex_descriptor place_type;
  typedef boost::graph_traits<PetriGraph>::vertex_descriptor transition_type;
  // This gspn expects transitions to be of this base class.
  // Derive from this base class to make transitions.
  typedef smv::ExplicitTransition<Local,RNG,ExtraState> Transition;

private:
  // The GSPN gets the marking type from the State.
  BiMap _bimap;
  std::map<transition_type,std::unique_ptr<Transition>> transitions;
  PetriGraph graph;

public:
  ExplicitTransitions(size_t num_vertices) : graph(num_vertices) {}

  // These copy constructor and copy assignment operators are
  // deleted because this type can only be moved with std::move();
  ExplicitTransitions(const ExplicitTransitions&)=delete;
  ExplicitTransitions& operator=(const ExplicitTransitions&)=delete;


  ExplicitTransitions(ExplicitTransitions&& other)
  {
    transitions=std::move(other.transitions);
    graph=std::move(other.graph);
    _bimap=std::move(other._bimap);
  }


  ExplicitTransitions& operator=(ExplicitTransitions&& other)
  {
    if (this!=other)
    {
      transitions=std::move(other.transitions);
      graph=std::move(other.graph);
      _bimap=std::move(other._bimap);
    }
  }

  ~ExplicitTransitions() {}

  size_t place_vertex(PlaceKey p) const
  {
    return get_pvertex(_bimap, p);
  }


  PlaceKey vertex_place(size_t v) const
  {
    return get_place(_bimap, v);
  }


  size_t transition_vertex(TransitionKey t) const
  {
    return get_tvertex(_bimap, t);
  }


  TransitionKey vertex_transition(size_t v) const
  {
    return get_transition(_bimap, v);
  }


  friend BuildGraph<ExplicitTransitions<PKey,TKey,Local,ETRand,ExtraState>>;

  template<typename State, typename P, typename T, typename L, typename Random>
  friend
  std::vector<std::tuple<size_t,size_t,int>>
  neighbors_of_transition(
    ExplicitTransitions<State,P,T,L,Random>& et,
    typename ExplicitTransitions<State,P,T,L,Random>::transition_type trans_id);

  template<typename F, typename State, typename P, typename T, typename L,
    typename Random>
  friend
  void neighbors_of_places(
    ExplicitTransitions<State,P,T,L,Random>& et,
    const std::set<typename ExplicitTransitions<State,P,T,L,Random>::place_type>&
    place_id, F func);

  template<typename Transitions, typename... Args>
  friend
  std::pair<bool,std::unique_ptr<TransitionDistribution<typename Transitions::RNG>>>
  enabled(const Transitions& et, typename Transitions::transition_type trans_id,
    Args&&... args);

  template<typename Transitions, typename... Args>
  friend
  void
  fire(Transitions& et, typename Transitions::transition_type trans_id,
    Args&&... args);

};



// Now that we have the object, this is how it fulfills the
// GSPN concept.
template<typename P, typename T, typename L, typename Random, typename State>
struct petri_place<ExplicitTransitions<P,T,L,Random,State>>
{
  typedef typename ExplicitTransitions<P,T,L,Random,State>::place_type type;
};



template<typename P, typename T, typename L, typename Random, typename State>
struct petri_transition<ExplicitTransitions<P,T,L,Random,State>>
{
  typedef typename ExplicitTransitions<P,T,L,Random,State>::transition_type type;
};



template<typename P, typename T, typename L, typename Random, typename State>
std::vector<std::tuple<size_t,size_t,int>>
neighbors_of_transition(
  ExplicitTransitions<P,T,L,Random,State>& et,
  typename ExplicitTransitions<P,T,L,Random,State>::transition_type trans_id)
{
  return neighbors_of_transition(et.graph, trans_id);
}


template<typename F, typename P, typename T, typename L,
  typename Random, typename State>
void neighbors_of_places(
  ExplicitTransitions<P,T,L,Random,State>& et,
  const std::set<typename ExplicitTransitions<P,T,L,Random,State>::place_type>&
  place_id, F func)
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
