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
#include "build_graph.hpp"
#include "gspn.hpp"
#include "logging.hpp"


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
  Enabled(const ExtraState& s, const LM& lm,
          double enabling_time, double current_time) const {
    BOOST_LOG_TRIVIAL(debug) << "The base enabled is unlikely correct to call";
    return {false,std::unique_ptr<TransitionDistribution<RNG>>(nullptr)};
  }

  virtual void Fire(ExtraState& s, LM& lm, RNG& rng) const {
    lm.template TransferByStochiometricCoefficient<0>(rng);
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
  typedef PKey UserPlaceKey;
  typedef TKey UserTransitionKey;
  typedef ETRand RNG;
  typedef boost::graph_traits<PetriGraph>::vertex_descriptor PlaceKey;
  typedef boost::graph_traits<PetriGraph>::vertex_descriptor TransitionKey;
  // This gspn expects transitions to be of this base class.
  // Derive from this base class to make transitions.
  typedef smv::ExplicitTransition<Local,RNG,ExtraState> Transition;

private:
  // The GSPN gets the marking type from the State.
  BiMap _bimap;
  std::map<TransitionKey,std::unique_ptr<Transition>> transitions;
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

  size_t PlaceVertex(UserPlaceKey p) const
  {
    return GetPvertex(_bimap, p);
  }


  UserPlaceKey VertexPlace(size_t v) const
  {
    return GetPlace(_bimap, v);
  }


  size_t TransitionVertex(UserTransitionKey t) const
  {
    return GetTvertex(_bimap, t);
  }


  UserTransitionKey VertexTransition(size_t v) const
  {
    return GetTransition(_bimap, v);
  }


  friend BuildGraph<ExplicitTransitions<PKey,TKey,Local,ETRand,ExtraState>>;

  template<typename State, typename P, typename T, typename L, typename Random>
  friend
  std::vector<std::tuple<size_t,size_t,int>>
  NeighborsOfTransition(
    ExplicitTransitions<State,P,T,L,Random>& et,
    typename ExplicitTransitions<State,P,T,L,Random>::TransitionKey trans_id);

  template<typename F, typename State, typename P, typename T, typename L,
    typename Random>
  friend
  void NeighborsOfPlaces(
    ExplicitTransitions<State,P,T,L,Random>& et,
    const std::set<typename ExplicitTransitions<State,P,T,L,Random>::PlaceKey>&
    place_id, F func);

  template<typename Transitions, typename... Args>
  friend
  std::pair<bool,std::unique_ptr<TransitionDistribution<typename Transitions::RNG>>>
  Enabled(const Transitions& et, typename Transitions::TransitionKey trans_id,
    Args&&... args);

  template<typename Transitions, typename... Args>
  friend
  void Fire(Transitions& et, typename Transitions::TransitionKey trans_id,
    Args&&... args);
};



template<typename P, typename T, typename L, typename Random, typename State>
std::vector<std::tuple<size_t,size_t,int>>
NeighborsOfTransition(
  ExplicitTransitions<P,T,L,Random,State>& et,
  typename ExplicitTransitions<P,T,L,Random,State>::TransitionKey trans_id)
{
  return NeighborsOfTransition(et.graph, trans_id);
}


template<typename F, typename P, typename T, typename L,
  typename Random, typename State>
void NeighborsOfPlaces(
  ExplicitTransitions<P,T,L,Random,State>& et,
  const std::set<typename ExplicitTransitions<P,T,L,Random,State>::PlaceKey>&
  place_id, F func)
{
  return NeighborsOfPlaces(et.graph, place_id, func);
}


template<typename Transitions, typename... Args>
std::pair<bool,std::unique_ptr<TransitionDistribution<typename Transitions::RNG>>>
Enabled(const Transitions& et, typename Transitions::TransitionKey trans_id,
  Args&&... args)
{
  return et.transitions.at(trans_id)->Enabled(std::forward<Args>(args)...);
}



template<typename Transitions, typename... Args>
void Fire(Transitions& et, typename Transitions::TransitionKey trans_id,
  Args&&... args)
{
  et.transitions.at(trans_id)->Fire(std::forward<Args>(args)...);
}

}
}
#endif // _EXPLICIT_TRANSITIONS_H_
