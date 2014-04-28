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
#include <tuple>
#include <map>
#include <iostream>
#include <limits>
#include <algorithm>
#include <memory>
#include <random>
#include <set>
#include <functional>
#include "stochnet.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/mpl/vector.hpp"
#include "smv.hpp"
#include "cow_token.hpp"


namespace smv=afidd::smv;
using namespace afidd::smv;

using CowGen=std::mt19937;




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
  CowPlace(int d, int i, int s) : disease(d), individual(i), subgroup(s) {}

  friend inline
  bool operator<(const CowPlace& a, const CowPlace& b) {
    return LazyLess(a.disease, b.disease, a.individual,
      b.individual, a.subgroup, b.subgroup);
  }


  friend inline
  bool operator==(const CowPlace& a, const CowPlace& b) {
    return (a.disease==b.disease)&& (a.individual==b.individual)
      && (a.subgroup==b.subgroup);
  }


  friend inline
  std::ostream&
  operator<<(std::ostream& os, const CowPlace& cp) {
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
  bool operator<(const CowT& a, const CowT& b) {
    return LazyLess(a.cow1, b.cow1, a.cow2, b.cow2,
      a.sg1, b.sg1, a.sg2, b.sg2, a.kind, b.kind);
  }

  friend inline
  bool operator==(const CowT& a, const CowT& b) {
    return (a.cow1==b.cow1) && (a.cow2==b.cow2)
        && (a.sg1==b.sg1) && (a.sg2==b.sg2) && (a.kind==b.kind);
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const CowT& cp) {
    return os << '(' << cp.cow1 << "," << cp.cow2 << ","
      << cp.sg1 << ',' << cp.sg2 << ',' << ", " << cp.kind << ')';
  }
};




// This is the central set of definitions in order to use the
// ExplicitTransitions representation of the GSPN.
using Local=smv::LocalMarking<smv::Colored<Cow>,
    smv::Uncolored<std::map<int,double>>>;
using CowTransitions=smv::ExplicitTransitions<
    CowPlace,CowT,Local,CowGen>;

// We make a transition that meets the requirements of the GSPN object
// by deriving it from the Transition type it defines.

class CowTransition
: public smv::ExplicitTransition<Local,CowGen>
{
public:
  CowTransition(int64_t cow_id) : cow_id(cow_id) {}
  virtual ~CowTransition() {}

  int64_t cow_id;
};


// These are just shorthand.
using Dist=smv::TransitionDistribution<CowGen>;
using ExpDist=smv::ExponentialDistribution<CowGen>;



class InfectNeighbor : public CowTransition
{
public:
  InfectNeighbor(int64_t cow_id) : CowTransition(cow_id) {}

  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0) const {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) const {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};







/*! The Petri Net we make depends only on the local marking, not
 *  the marking, because transitions are defined on local state.
 */
CowTransitions
Herd(int64_t initial_cnt, int64_t total_cnt)
{
  BuildGraph<CowTransitions> bg;
  using PlaceEdge=BuildGraph<CowTransitions>::PlaceEdge;
  // subgroups
  enum { c, h1, h2, d, death, sale, culling};
  // disease states
  enum { dm, ds, dti, dr, dpi};

  for (auto sg : {c, h1, h2, d, death, sale, culling}) {
    for (int who=0; who<total_cnt; ++who) {
      for (auto disease : std::vector<int>{dm, ds, dti, dr, dpi}) {
        bg.AddPlace({disease, who, sg}, 0);
      }
      // same group, same cow, kind=0 is becoming susceptible.
      bg.AddTransition({who, who, sg, sg, 0},
        {PlaceEdge{CowPlace{dm, who, sg}, -1},
         PlaceEdge{CowPlace{ds, who, sg}, 1}},
        std::unique_ptr<CowTransition>(new InfectNeighbor(who))
        );

      // same group, same cow, kind=1 is recovering.
      bg.AddTransition({who, who, sg, sg, 1},
        {PlaceEdge{CowPlace{dti, who, sg}, -1},
         PlaceEdge{CowPlace{dr, who, sg}, 1}},
        std::unique_ptr<CowTransition>(new InfectNeighbor(who))
        );
    }
  }

  for (auto sg : {c, h1, h2, d, death, sale, culling}) {
    for (int who=0; who<total_cnt; ++who) {
      for (auto disease : std::vector<int>{dm, ds, dti, dr, dpi})
      {
      }
    }
  }
  
  return std::move(bg.Build());
}




int main(int argc, char *argv[])
{
  CowGen rng{1};

  auto gspn=Herd(100, 10);

  using Mark=smv::Marking<int64_t,
      smv::Colored<Cow>,smv::Uncolored<std::map<int,double>>>;
  Mark m;
  using CowState=smv::GSPNState<Mark,int64_t>;
  CowState state;

  enum { lambda, beta, gamma };
  std::map<int,double> params;
  params[lambda]=1.0;
  params[beta]=1.7;
  params[gamma]=0.5;
  int64_t params_place_id=0;
  Add<1>(m, params_place_id, params);
  assert(Length<1>(m, params_place_id)==1);
  assert(Length<0>(m, 27, 13)==0);

  using Propagator=PropagateCompetingProcesses<int64_t,CowGen>;
  using Dynamics=StochasticDynamics<CowTransitions,CowState,CowGen>;
  Propagator competing;
  Dynamics dynamics(gspn, {&competing});
  dynamics.Initialize(&state, &rng);
  bool running=dynamics(state);
  return 0;
}

