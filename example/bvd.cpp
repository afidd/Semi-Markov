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
#include "bvd_distributions.hpp"
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
  int kind;

  CowT()=default;
  CowT(int k): kind(k){}

  friend inline
  bool operator<(const CowT& a, const CowT& b) {
    return  a.kind < b.kind;
  }

  friend inline
  bool operator==(const CowT& a, const CowT& b) {
    return a.kind==b.kind;
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const CowT& cp) {
    return os << '(' << cp.kind << ')';
  }
};


struct CowUserState {
  std::map<std::string,double> params;
};



// This is the central set of definitions in order to use the
// ExplicitTransitions representation of the GSPN.
using Local=smv::LocalMarking<smv::Colored<Cow>,
    smv::Uncolored<std::map<int,double>>>;
using CowTransitions=smv::ExplicitTransitions<
    CowPlace,CowT,Local,CowGen,CowUserState>;

// We make a transition that meets the requirements of the GSPN object
// by deriving it from the Transition type it defines.

using CowTransition=smv::ExplicitTransition<Local,CowGen,CowUserState>;

// These are just shorthand.
using Dist=smv::TransitionDistribution<CowGen>;
using ExpDist=smv::ExponentialDistribution<CowGen>;




class CalfToH1 : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      double weaning=s.params.at("weaning");
      auto dist=new smv::DiracDistribution<CowGen>(weaning, te);
      return {true, std::unique_ptr<Dist>(dist)};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};



class H1ToH2 : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      double left=s.params.at("breeding_left");
      double middle=s.params.at("breeding_middle");
      double right=s.params.at("breeding_right");
      return {true, std::unique_ptr<Dist>(
          new smv::TriangularDistribution<CowGen>(left, middle, right, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};



class SellCalf : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};



class H2ToDairy : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<Dist>(
        new HeiferCalvingDistribution<CowGen>(te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};



class Parturition : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};



class CalfToDeath : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};



class H2ToDeath : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};



class DairyToDeath : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};




class H2Cull : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};




class DairyCull : public CowTransition
{
public:
  virtual std::pair<bool,std::unique_ptr<TransitionDistribution<CowGen>>>
  Enabled(const UserState& s, const Local& lm, double te, double t0, CowGen& rng) {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0, CowGen& rng) {
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
  enum { c, h1, h2, d, death, sold, culling};
  // disease states
  enum { dm, ds, dti, dr, dpi};
  // Transition kinds
  enum { moveherd, purchase, sale, cull, die, birth };

  for (auto sg : std::vector<int>{c, h1, h2, d, death, sale, culling}) {
    for (int who=0; who<total_cnt; ++who) {
      for (auto disease : std::vector<int>{dm, ds, dti, dr, dpi}) {
        bg.AddPlace(CowPlace{disease, who, sg}, 0);
      }
    }
  }

  for (int who=0; who<total_cnt; ++who) {
    bg.AddTransition({moveherd},
      {PlaceEdge{CowPlace{ds, who, c}, -1},
       PlaceEdge{CowPlace{ds, who, h1}, 1}},
       std::unique_ptr<CowTransition>(new CalfToH1()));
  }
  for (int who=0; who<total_cnt; ++who) {
    bg.AddTransition({moveherd},
      {PlaceEdge{CowPlace{ds, who, h1}, -1},
       PlaceEdge{CowPlace{ds, who, h2}, 1}},
       std::unique_ptr<CowTransition>(new H1ToH2()));
  }
  for (int who=0; who<total_cnt; ++who) {
    bg.AddTransition({sale},
      {PlaceEdge{CowPlace{ds, who, c}, -1},
       PlaceEdge{CowPlace{ds, who, sold}, 1}},
       std::unique_ptr<CowTransition>(new SellCalf()));
  }
  // Why? A new calf's identity depends on which slot is free.
  // Cannot determine it at time of GSPN construction.
  std::vector<PlaceEdge> anycalf(total_cnt+2);
  for (int destination=0; destination<total_cnt; ++destination) {
    anycalf[2+destination]=PlaceEdge{CowPlace{ds, destination, c}, 1};
  }
  for (int who=0; who<total_cnt; ++who) {
    anycalf[0]=PlaceEdge{CowPlace{ds, who, h2}, -1};
    anycalf[1]=PlaceEdge{CowPlace{ds, who, d}, 1};
    bg.AddTransition({birth}, anycalf,
       std::unique_ptr<CowTransition>(new H2ToDairy()));
  }
  for (int who=0; who<total_cnt; ++who) {
    anycalf[0]=PlaceEdge{CowPlace{ds, who, d}, -1};
    anycalf[1]=PlaceEdge{CowPlace{ds, who, d}, 1};
    bg.AddTransition({birth}, anycalf,
       std::unique_ptr<CowTransition>(new Parturition()));
  }

  // perinatal death
  for (int who=0; who<total_cnt; ++who) {
    bg.AddTransition({die},
      {PlaceEdge{CowPlace{ds, who, c}, -1},
       PlaceEdge{CowPlace{ds, who, death}, 1}},
       std::unique_ptr<CowTransition>(new CalfToDeath()));
  }

  // no H1 death, but H2 death
  for (int who=0; who<total_cnt; ++who) {
    bg.AddTransition({die},
      {PlaceEdge{CowPlace{ds, who, h2}, -1},
       PlaceEdge{CowPlace{ds, who, death}, 1}},
       std::unique_ptr<CowTransition>(new H2ToDeath()));
  }
  
  // culling of dairy
  for (int who=0; who<total_cnt; ++who) {
    bg.AddTransition({die},
      {PlaceEdge{CowPlace{ds, who, d}, -1},
       PlaceEdge{CowPlace{ds, who, death}, 1}},
       std::unique_ptr<CowTransition>(new DairyToDeath()));
  }
  

  for (int who=0; who<total_cnt; ++who) {
    bg.AddTransition({cull},
      {PlaceEdge{CowPlace{ds, who, h2}, -1},
       PlaceEdge{CowPlace{ds, who, culling}, 1}},
       std::unique_ptr<CowTransition>(new H2Cull()));
  }
  for (int who=0; who<total_cnt; ++who) {
    bg.AddTransition({cull},
      {PlaceEdge{CowPlace{ds, who, d}, -1},
       PlaceEdge{CowPlace{ds, who, culling}, 1}},
       std::unique_ptr<CowTransition>(new DairyCull()));
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
  using CowState=smv::GSPNState<Mark,int64_t,CowUserState>;
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

