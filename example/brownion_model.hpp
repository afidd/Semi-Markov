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
#ifndef _BROWNION_MODEL_H_
#define _BROWNION_MODEL_H_ 1
#include <random>
#include "gspn.hpp"
#include "continuous_state.hpp"
#include "marking.hpp"
#include "distributions.hpp"
#include "smv_algorithm.hpp"
#include "local_marking.hpp"


namespace smv=afidd::smv;
using namespace smv;

using RandGen=std::mt19937;





struct PlaceType
{
  PlaceType()=default;
  PlaceType(int64_t i, int64_t j, int state) : i(i), j(j), state(state) {}
  PlaceType(const PlaceType&)=default;
  PlaceType& operator=(const PlaceType&)=default;

  int64_t i, j;
  int state;

  friend inline
  std::ostream& operator<<(std::ostream& os, const PlaceType& id) {
    return os << "{" << id.i << "," << id.j << "," << id.state << ")";
  }

  friend inline
  bool operator==(const PlaceType& a, const PlaceType& b) {
    return (a.i==b.i) && (a.j==b.j) && (a.state==b.state);
  }

  friend inline
  bool operator<(const PlaceType& a, const PlaceType& b) {
    return LazyLess(a.i, b.i, a.j, b.j, a.state, b.state);
  }
  
};



struct TransitionType
{
  TransitionType()=default;
  TransitionType(PlaceType from, PlaceType to) : from(from), to(to) {}
  TransitionType(const TransitionType&)=default;
  TransitionType& operator=(const TransitionType&)=default;
  PlaceType from, to;

  friend inline
  std::ostream& operator<<(std::ostream& os, const TransitionType& id) {
    return os << "(" << id.from << "," << id.to << ")";
  }

  friend inline
  bool operator<(const TransitionType& a, const TransitionType& b) {
    return LazyLess(a.from, b.from, a.to, b.to);
  }

  friend inline
  bool operator==(const TransitionType& a, const TransitionType& b) {
    return (a.from==b.from) && (a.to==b.to);
  }
};


struct UserState {};


struct IndividualToken {};

using Local=smv::LocalMarking<smv::Uncolored<IndividualToken>>;
using Mark=smv::Marking<PlaceType,
                        smv::Uncolored<IndividualToken>>;
using BrownionState=smv::GSPNState<Mark,UserState>;

using Dist=smv::TransitionDistribution<RandGen>;
using ExpDist=smv::ExponentialDistribution<RandGen>;
using Weibull=smv::WeibullDistribution<RandGen>;


class BrownionGSPN
{
public:
  typedef PlaceType PlaceKey;
  typedef TransitionType TransitionKey;
  // Could store the state parameters and distributions here
  // if we wanted.
};



std::pair<bool,std::unique_ptr<TransitionDistribution<RandGen>>>
Enabled(const BrownionGSPN& et, TransitionType trans_id,
  const UserState& s, const Local& lm, double te, double t0)
{
  if (lm.template Length<0>(0)>0) {
    // This is where we choose the distributions for the two
    // Brownion states.
    if (trans_id.from.state==0) {
      return {true, std::unique_ptr<Weibull>(new Weibull(1.0,1.0, te))};
    } else {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    }
  } else {
    return {false, std::unique_ptr<Dist>(nullptr)};
  }
}




template<typename RNG>
void
Fire(BrownionGSPN& et, TransitionType trans_id,
  UserState& s, Local& lm, RNG& rng)
{
  lm.template Move<0,0>(0, 1, 1);
}



std::vector<std::tuple<PlaceType,int,int>>
NeighborsOfTransition(BrownionGSPN& g, TransitionType trans_id)
{
  std::vector<std::tuple<PlaceType,int,int>> place_ids;
  place_ids.push_back(std::make_tuple(trans_id.from, 0, -1));
  place_ids.push_back(std::make_tuple(trans_id.to, 0, -1));
  return place_ids;
}



template<typename F>
void NeighborsOfPlaces(BrownionGSPN& g,
  const std::set<PlaceType>& place_id, const F& func)
{
  for (auto p : place_id) {
    // Transitions that start at this place.
    func(TransitionType{p, {p.i, p.j-1, p.state}});
    func(TransitionType{p, {p.i, p.j+1, p.state}});
    func(TransitionType{p, {p.i-1, p.j, p.state}});
    func(TransitionType{p, {p.i+1, p.j, p.state}});

    // Transitions that end at this place.
    func(TransitionType{{p.i+1, p.j, p.state}, p});
    func(TransitionType{{p.i-1, p.j, p.state}, p});
    func(TransitionType{{p.i, p.j+1, p.state}, p});
    func(TransitionType{{p.i, p.j-1, p.state}, p});

    // Transition to change Brownion state.
    func(TransitionType{p, {p.i, p.j, 1-p.state}});
  }
}


#endif /* _BROWNION_MODEL_H_ */
