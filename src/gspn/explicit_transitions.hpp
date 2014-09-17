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
#include "distributions.hpp"
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
          double enabling_time, double current_time) {
    SMVLOG(BOOST_LOG_TRIVIAL(debug) <<
      "The base enabled is unlikely correct to call");
    return {false,std::unique_ptr<TransitionDistribution<RNG>>(nullptr)};
  }

  virtual void Fire(ExtraState& s, LM& lm, double t0, RNG& rng) {
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
  typedef BiGraphCorrespondence<PKey,TKey,int64_t> BiMap;

public:
  typedef PKey UserPlaceKey;
  typedef TKey UserTransitionKey;
  typedef ETRand RNG;
  typedef int64_t PlaceKey;
  typedef int64_t TransitionKey;
  // This gspn expects transitions to be of this base class.
  // Derive from this base class to make transitions.
  typedef smv::ExplicitTransition<Local,RNG,ExtraState> Transition;

private:
  // The GSPN gets the marking type from the State.
  BiMap bimap_;
  // Transitions are stored as pointers on the free store because they
  // are polymorphic and base class information would be lost if they
  // were stored as bare objects.
  std::map<TransitionKey,std::unique_ptr<Transition>> transitions;
  PetriGraph graph;

public:
  explicit ExplicitTransitions(size_t num_vertices) : graph(num_vertices) {}

  // These copy constructor and copy assignment operators are
  // deleted because this type can only be moved with std::move();
  ExplicitTransitions(const ExplicitTransitions&)=delete;
  ExplicitTransitions& operator=(const ExplicitTransitions&)=delete;


  ExplicitTransitions(ExplicitTransitions&& other)
  {
    transitions=std::move(other.transitions);
    graph=std::move(other.graph);
    bimap_=std::move(other.bimap_);
  }


  ExplicitTransitions& operator=(ExplicitTransitions&& other)
  {
    if (this!=other)
    {
      transitions=std::move(other.transitions);
      graph=std::move(other.graph);
      bimap_=std::move(other.bimap_);
    }
  }

  ~ExplicitTransitions() {}

  int64_t PlaceVertex(UserPlaceKey p) const
  {
    return GetPvertex(bimap_, p);
  }


  UserPlaceKey VertexPlace(int64_t v) const
  {
    return GetPlace(bimap_, v);
  }


  int64_t TransitionVertex(UserTransitionKey t) const
  {
    return GetTvertex(bimap_, t);
  }


  UserTransitionKey VertexTransition(int64_t v) const
  {
    return GetTransition(bimap_, v);
  }

  // Given a transition, find all the places that neighbor it.
  // This exists so that we can ask questions about a transition that fired.
  std::tuple<int64_t,UserPlaceKey>
  PlaceOfTransition(int64_t transition_id, int examine_edge) const {
    int edge_idx=0;
    auto initer=in_edges(transition_id, graph);
    for ( ; initer.first!=initer.second; ++initer.first) {
      if (edge_idx==examine_edge) {
        int64_t place_vertex=source(*initer.first, graph);
        auto place_key=GetPlace(bimap_, place_vertex);
        return std::make_tuple(place_vertex, place_key);
      }
      ++edge_idx;
    }
    auto outiter=out_edges(transition_id, graph);
    for (; outiter.first!=outiter.second; ++outiter.first) {
      if (edge_idx==examine_edge) {
        int64_t place_vertex=target(*outiter.first, graph);
        auto place_key=GetPlace(bimap_, place_vertex);
        return std::make_tuple(place_vertex, place_key);
      }
      ++edge_idx;
    }
    SMVLOG(BOOST_LOG_TRIVIAL(error)<<
      "Could not find a place near a transition.");
    assert(false);
    return std::make_tuple(0, UserPlaceKey{});
  }


  template<typename Writer>
  void WriteKeys(const Writer& writer) {
    for (const auto& kv : transitions) {
      writer.Write(kv.first, this->VertexTransition(kv.first));
    }
  }

  friend BuildGraph<ExplicitTransitions<PKey,TKey,Local,ETRand,ExtraState>>;

  template<typename State, typename P, typename T, typename L, typename Random>
  friend
  std::vector<std::tuple<int64_t,int,int>>
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
std::vector<std::tuple<int64_t,int,int>>
NeighborsOfTransition(
  ExplicitTransitions<P,T,L,Random,State>& et,
  typename ExplicitTransitions<P,T,L,Random,State>::TransitionKey trans_id)
{
  return GraphNeighborsOfTransition<PetriGraphType>(et.graph, trans_id);
}


template<typename F, typename P, typename T, typename L,
  typename Random, typename State>
void NeighborsOfPlaces(
  ExplicitTransitions<P,T,L,Random,State>& et,
  const std::set<typename ExplicitTransitions<P,T,L,Random,State>::PlaceKey>&
  place_id, F func)
{
  return GraphNeighborsOfPlaces(et.graph, place_id, func);
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
