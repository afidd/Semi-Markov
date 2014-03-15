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
#include "cow_token.h"


namespace smv=afidd::smv;
using namespace afidd::smv;

using CowGen=boost::random::mt19937;




using PN=smv::PetriGraphType;

namespace afidd{
namespace smv{
template<>
struct petri_place<PetriGraphType>
{
  typedef size_t type;
};


template<>
struct petri_transition<PetriGraphType>
{
  typedef size_t type;
};


template<>
struct petri_graph<PetriGraphType>
{
  typedef PetriGraphType type;
};
} // smv
} // afidd

using Mark=smv::Marking<smv::place_t<PN>,
    smv::Colored<Cow>,smv::Uncolored<std::map<size_t,double>>>;
using CowState=smv::GSPNState<PN,Mark>;

class CowTransition
: public smv::ExplicitTransition<smv::LocalMarking<Mark>,CowState,CowGen>
{
public:
  CowTransition(size_t cow_id) : cow_id(cow_id) {}
  virtual ~CowTransition() {}

  size_t cow_id;
};


using Dist=smv::TransitionDistribution<CowGen>;
using ExpDist=smv::ExponentialDistribution<CowGen>;
using NoDist=smv::NoDistribution<CowGen>;
using CowTransitions=smv::ExplicitTransitions<smv::LocalMarking<Mark>,
    CowState,PN,CowGen>;

class InfectNeighbor : public CowTransition
{
public:
  InfectNeighbor(size_t cow_id) : CowTransition(cow_id) {}

  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  enabled(const CowState& s,
    const LocalMarking<Mark>& lm) const
  {
    return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0))};
  }

  virtual void fire(CowState& s, LocalMarking<Mark>& lm, CowGen& rng) const
  {
    return;
  }
};



/*! The Petri Net we make depends only on the local marking, not
 *  the marking, because transitions are defined on local state.
 */
template<typename PN, typename LocalMarking, typename CowState>
CowTransitions
herd(size_t initial_cnt, size_t total_cnt)
{
  enum { c, h1, h2, d, death, sale, culling};
  graph_t<PN> g(7);
  std::map<size_t,double> params;

  PetriGraphVertexProperty vprop;
  for (auto place : std::vector<int>{c, h1, h2, d, death, sale, culling})
  {
    vprop.color=PetriGraphColor::Place;
    vprop.token_layer=0;
    g[place]=vprop;
  }

  vprop.color=PetriGraphColor::Transition;

  auto et=CowTransitions(g);
  et.transitions.emplace(17, std::move(std::unique_ptr<CowTransition>(
      new InfectNeighbor(17))));

  return std::move(et);
}




int main(int argc, char *argv[])
{
  Mark m;
  enum { lambda, beta, gamma };
  std::map<size_t,double> params;
  params[lambda]=1.0;
  params[beta]=1.7;
  params[gamma]=0.5;
  place_t<PN> params_place_id=0;
  add<1>(m, params_place_id, params);
  assert(length<1>(m, params_place_id)==1);
  assert(length<0>(m, 27, 13)==0);

  CowGen rng{1};

  CowState state;
  auto gspn=herd<PN,LocalMarking<Mark>,CowState>(100, 10);

  PartialCoreMatrix<CowTransitions,CowState,CowGen>
      system(gspn, state);
  auto token=[](CowState&) { };
  auto next=propagate_competing_processes(system, token, rng);
  if (std::get<1>(next)<=std::numeric_limits<double>::max())
  {
    std::cout << "We have a value." << std::endl;
  }
  return 0;
}

