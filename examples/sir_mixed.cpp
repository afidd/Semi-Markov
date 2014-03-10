/*! Fully mixed SIR with arbitrary transition rates.
 */

#include <tuple>
#include <map>
#include <iostream>
#include <limits>
#include <algorithm>
#include <memory>
#include <set>
#include <functional>
#include "stochnet.h"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/mpl/vector.hpp"
#include "gspn.h"
#include "petri_graph.h"
#include "marking.h"
#include "distributions.h"
#include "continuous_state.h"
#include "explicit_transitions.h"
#include "partial_core_matrix.h"
#include "continuous_dynamics.h"
#include "logging.h"

namespace smv=afidd::smv;
using namespace smv;
using RandGen=boost::random::mt19937;




using PG=PetriGraphType;


// Inject traits types for into the afidd namespace where
// they are used to determine place and transition token types.
namespace afidd{
namespace smv{
template<>
struct petri_place<PG>
{
  typedef size_t type;
};


template<>
struct petri_transition<PG>
{
  typedef size_t type;
};


template<>
struct petri_graph<PG>
{
  typedef PG type;
};
} // smv
} // afidd


struct IndividualToken
{
  IndividualToken()=default;
};


// Marking of the net.
using Mark=Marking<place_t<PG>, Uncolored<IndividualToken>>;
// State of the continuous dynamical system.
using SIRState=GSPNState<PG, Mark>;

class SIRTransition
: public ExplicitTransition<LocalMarking<Mark>, SIRState, RandGen>
{
public:
  SIRTransition() {}
  virtual ~SIRTransition() {}
};

// This class holds the transitions.
using SIRGSPN=
    ExplicitTransitions<LocalMarking<Mark>,SIRState, PG, RandGen>;


using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;
using NoDist=NoDistribution<RandGen>;




// Now make specific transitions.
class InfectNeighbor : public SIRTransition
{

  virtual std::pair<bool, std::unique_ptr<Dist>>
  enabled(const SIRState& s, const LocalMarking<Mark>& lm) const override
  {
    if (lm.template input_tokens_sufficient<0>())
    {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0))};
    }
    else
    {
      return {false, std::unique_ptr<NoDist>(new NoDist())};
    }
  }

  virtual void fire(SIRState& s, LocalMarking<Mark>& lm,
      RandGen& rng) const override
  {
    BOOST_LOG_TRIVIAL(debug) << "Fire infection " << lm;
    lm.template transfer_by_stochiometric_coefficient<0>(rng);
  }

};





// Now make specific transitions.
class Recover : public SIRTransition
{

  virtual std::pair<bool, std::unique_ptr<Dist>>
  enabled(const SIRState& s, const LocalMarking<Mark>& lm) const override
  {
    if (lm.template input_tokens_sufficient<0>())
    {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0))};
    }
    else
    {
      return {false, std::unique_ptr<NoDist>(new NoDist())};
    }
  }

  virtual void fire(SIRState& s, LocalMarking<Mark>& lm,
      RandGen& rng) const override
  {
    BOOST_LOG_TRIVIAL(debug) << "Fire recovery "<< lm;
    lm.template transfer_by_stochiometric_coefficient<0>(rng);
  }

};



/*! SIR infection on an all-to-all graph of uncolored tokens.
 */
SIRGSPN
build_system(size_t individual_cnt)
{
  size_t individual_state_cnt=3;
  size_t place_cnt=individual_cnt * individual_state_cnt;
  // Two infections for each combination, so the twos cancel.
  size_t infection_cnt=individual_cnt*(individual_cnt-1);
  size_t recovery_cnt=individual_cnt;
  size_t transition_cnt=infection_cnt+recovery_cnt;

  PG build_graph(place_cnt + transition_cnt);
  SIRGSPN et(build_graph);
  PG& graph=et.graph;

  enum { s, i, r };
  PetriGraphVertexProperty vprop;
  size_t assigned_place_idx;
  for (auto ind_idx=0; ind_idx<individual_cnt; ind_idx++)
  {
    for (auto place : std::vector<int>{s, i, r})
    {
      vprop.color=PetriGraphColor::Place;
      vprop.token_layer=0;
      assigned_place_idx=ind_idx*individual_state_cnt+place;
      graph[assigned_place_idx]=vprop;
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "Places from "<<0<<" to "<<place_cnt;
  assert(assigned_place_idx=place_cnt-1);

  vprop.color=PetriGraphColor::Transition;

  size_t trans_idx=place_cnt;
  for (auto left_idx=0; left_idx<individual_cnt-1; left_idx++)
  {
    for (auto right_idx=left_idx+1; right_idx<individual_cnt; right_idx++)
    {
      graph[trans_idx]=vprop;
      add_edge(left_idx*3+i, trans_idx, {-1}, graph);
      add_edge(right_idx*3+s, trans_idx, {-1}, graph);
      add_edge(trans_idx, left_idx*3+i, {1}, graph);
      add_edge(trans_idx, right_idx*3+i, {1}, graph);
      et.transitions.emplace(trans_idx++,
        std::move(std::unique_ptr<SIRTransition>(new InfectNeighbor())));

      graph[trans_idx]=vprop;
      add_edge(left_idx*3+s, trans_idx, {-1}, graph);
      add_edge(right_idx*3+i, trans_idx, {-1}, graph);
      add_edge(trans_idx, left_idx*3+i, {1}, graph);
      add_edge(trans_idx, right_idx*3+i, {1}, graph);
      et.transitions.emplace(trans_idx++,
        std::move(std::unique_ptr<SIRTransition>(new InfectNeighbor())));
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "Infections from "<<place_cnt<<" to "<<trans_idx;
  assert(trans_idx==place_cnt+infection_cnt);

  for (auto rec_idx=0; rec_idx<individual_cnt; rec_idx++)
  {
    graph[trans_idx]=vprop;
    add_edge(rec_idx*3+i, trans_idx, {-1}, graph);
    add_edge(trans_idx, rec_idx*3+r, {1}, graph);
      et.transitions.emplace(trans_idx++,
        std::move(std::unique_ptr<SIRTransition>(new Recover())));
  }
  BOOST_LOG_TRIVIAL(trace) << "Last transition "<< trans_idx-1;

  // std::move the transitions because they contain unique_ptr.
  return std::move(et);
}



std::ostream& operator<<(std::ostream& os, const Mark& m)
{
  const auto& mmap=std::get<0>(m._maps);

  for (auto kv : mmap)
  {
    os << "("<< kv.first << "," << kv.second.size() <<") ";
  }
  return os;
}




int main(int argc, char *argv[])
{
  afidd::log_init("debug");

  size_t individual_cnt=100;
  enum { beta, gamma };
  std::map<size_t,double> params;
  params[beta]=1.0;
  params[gamma]=0.5;

  RandGen rng(1);

  auto gspn=build_system(individual_cnt);

  SIRState state;
  for (size_t individual=0; individual<individual_cnt; ++individual)
  {
    add<0>(state.marking, 3*individual, IndividualToken{});
  }

  using Markov=PartialCoreMatrix<SIRGSPN, SIRState, RandGen>;
  Markov system(gspn, state);

  BOOST_LOG_TRIVIAL(debug) << state.marking;

  // The initial input string moves a token from susceptible to infected.
  auto first_case=smv::uniform_index(rng, individual_cnt);
  auto input_string=[&first_case](SIRState& state)->void {
    move<0,0>(state.marking, first_case*3, first_case*3+1, 1);
  };
  auto next=propagate_competing_processes(system, input_string, rng);

  auto nothing=[](SIRState&)->void {};
  for ( ;
    std::get<1>(next)<std::numeric_limits<double>::max();
    next=propagate_competing_processes(system, nothing, rng))
  {
    BOOST_LOG_TRIVIAL(debug) << "trans " << std::get<0>(next) << " time " <<
        std::get<1>(next);
    BOOST_LOG_TRIVIAL(debug) << state.marking;
  }

  return 0;
}

