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
#include "build_graph.h"
#include "smv_algorithm.h"


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
    const LocalMarking<Mark>& lm, double current_time) const
  {
    return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, current_time))};
  }

  virtual void fire(CowState& s, LocalMarking<Mark>& lm, CowGen& rng) const
  {
    return;
  }
};




/*! A POD type for the id to a cow place.
 *  This is what it takes to make a class that is simple but
 *  also can be a key in a dictionary and initialize with
 *  brace initializers.
 */
struct CowPlace
{
  int disease;
  int individual;
  int subgroup;

  CowPlace()=default;
  CowPlace(int d, int i, int s)
  : disease(d), individual(i), subgroup(s)
  {}

  friend inline
  bool operator<(const CowPlace& a, const CowPlace& b)
  {
    return lazy_less(a.disease, b.disease, a.individual,
      b.individual, a.subgroup, b.subgroup);
  }


  friend inline
  bool operator==(const CowPlace& a, const CowPlace& b)
  {
    return (a.disease==b.disease)&& (a.individual==b.individual)
      && (a.subgroup==b.subgroup);
  }


  friend inline
  std::ostream&
  operator<<(std::ostream& os, const CowPlace& cp)
  {
    return os << '(' << cp.disease << ", " << cp.individual
        << ", " << cp.subgroup << ')';
  }
};



/*! This identifies a cow transition.
 *  We are being luxurious. If it's a cow-to-cow infection,
 *  both cows identify the transition. Subgroup-to-subgroup
 *  can also be recorded. The kind is then an identifier for
 *  a particular infection or movement.
 */
struct CowT
{
  int cow1;
  int cow2;
  int sg1;
  int sg2;
  int kind;

  CowT()=default;
  CowT(int c1, int c2, int s1, int s2, int k)
  : cow1(c1), cow2(c2), sg1(s1), sg2(s2), kind(k)
  {}

  friend inline
  bool operator<(const CowT& a, const CowT& b)
  {
    return lazy_less(a.cow1, b.cow1, a.cow2, b.cow2,
      a.sg1, b.sg1, a.sg2, b.sg2, a.kind, b.kind);
  }

  friend inline
  bool operator==(const CowT& a, const CowT& b)
  {
    return (a.cow1==b.cow1) && (a.cow2==b.cow2)
        && (a.sg1==b.sg1) && (a.sg2==b.sg2) && (a.kind==b.kind);
  }

  friend inline
  std::ostream&
  operator<<(std::ostream& os, const CowT& cp)
  {
    return os << '(' << cp.cow1 << "," << cp.cow2 << ","
      << cp.sg1 << ',' << cp.sg2 << ',' << ", " << cp.kind << ')';
  }
};





/*! The Petri Net we make depends only on the local marking, not
 *  the marking, because transitions are defined on local state.
 */
template<typename PN, typename LocalMarking, typename CowState>
CowTransitions
herd(size_t initial_cnt, size_t total_cnt)
{
  BuildGraph<CowPlace,CowT> bg;
  using PlaceEdge=BuildGraph<CowPlace,CowT>::PlaceEdge;
  // subgroups
  enum { c, h1, h2, d, death, sale, culling};
  // disease states
  enum { dm, ds, dti, dr, dpi};
  std::map<size_t,double> params;

  PetriGraphVertexProperty vprop;
  for (auto sg : {c, h1, h2, d, death, sale, culling})
  {
    for (int who=0; who<total_cnt; ++who)
    {
      for (auto disease : std::vector<int>{dm, ds, dti, dr, dpi})
      {
        bg.add_place({disease, who, sg}, 0);
      }
      bg.add_transition({who, who, sg, sg, 0},
        {PlaceEdge{CowPlace{dm, who, sg}, -1},
         PlaceEdge{CowPlace{ds, who, sg}, 1}});
    }
  }

  PetriGraphType g;
  BuildGraph<CowPlace,CowT>::BiMap b;
  std::tie(g, b)=bg.compile();

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

