#ifndef _BROWNION_MODEL_H_
#define _BROWNION_MODEL_H_ 1
#include "gspn.h"
#include "continuous_state.h"
#include "marking.h"
#include "distributions.h"


namespace smv=afidd::smv;
using namespace smv;

using RandGen=boost::random::mt19937;


// This helps make a correct less than operator for places and transitions.
// I've gotten this logic wrong too many times before out of laziness. Hence.
template<typename T>
bool lazy_less(const T& a, const T& b)
{
  return (a<b);
}


template<typename T, typename...Ts>
bool lazy_less(const T& a, const T& b, const Ts&... args)
{
  if (a<b)
  {
    return true;
  }
  else if (a==b)
  {
    return lazy_less(args...);
  }
  else
  {
    return false;
  }
}



struct PlaceType
{
  PlaceType()=default;
  PlaceType(int i, int j, int state) : i(i), j(j), state(state) {}
  PlaceType(const PlaceType&)=default;
  PlaceType& operator=(const PlaceType&)=default;

  int i, j;
  int state;

  friend inline
  std::ostream& operator<<(std::ostream& os, const PlaceType& id)
  {
    return os << "{" << id.i << "," << id.j << "," << id.state << ")";
  }

  friend inline
  bool operator==(const PlaceType& a, const PlaceType& b)
  {
    return (a.i==b.i) && (a.j==b.j) && (a.state==b.state);
  }

  friend inline
  bool operator<(const PlaceType& a, const PlaceType& b)
  {
    return lazy_less(a.i, b.i, a.j, b.j, a.state, b.state);
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
  std::ostream& operator<<(std::ostream& os, const TransitionType& id)
  {
    return os << "(" << id.from << "," << id.to << ")";
  }

  friend inline
  bool operator<(const TransitionType& a, const TransitionType& b)
  {
    return lazy_less(a.from, b.from, a.to, b.to);
  }

  friend inline
  bool operator==(const TransitionType& a, const TransitionType& b)
  {
    return (a.from==b.from) && (a.to==b.to);
  }
};



struct BrownionGraph
{
};


namespace afidd
{
  namespace smv
  {
  template<>
  struct petri_place<BrownionGraph>
  {
    typedef PlaceType type;
  };

  template<>
  struct petri_transition<BrownionGraph>
  {
    typedef TransitionType type;
  };
  }
}

struct IndividualToken {};


using Mark=smv::Marking<smv::place_t<BrownionGraph>, smv::Uncolored<IndividualToken>>;
using BrownionState=smv::GSPNState<BrownionGraph,Mark>;

using Dist=smv::TransitionDistribution<RandGen>;
using ExpDist=smv::ExponentialDistribution<RandGen>;
using Weibull=smv::WeibullDistribution<RandGen>;
using NoDist=smv::NoDistribution<RandGen>;


class BrownionGSPN
{
  // Could store the state parameters and distributions here
  // if we wanted.
};


namespace afidd
{
  namespace smv
  {
  template<>
  struct petri_place<BrownionGSPN>
  {
    typedef PlaceType type;
  };

  template<>
  struct petri_transition<BrownionGSPN>
  {
    typedef TransitionType type;
  };
  }
}


namespace afidd
{
namespace smv
{
std::pair<bool,std::unique_ptr<TransitionDistribution<RandGen>>>
enabled(const BrownionGSPN& et, TransitionType trans_id,
  const BrownionState& s, const smv::LocalMarking<Mark>& lm)
{
  if (lm.template length<0>(0)>0)
  {
    // This is where we choose the distributions for the two
    // Brownion states.
    if (trans_id.from.state==0)
    {
      return {true, std::unique_ptr<Weibull>(new Weibull(1.0,1.0))};
    }
    else
    {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0))};
    }
  }
  else
  {
    return {false, std::unique_ptr<NoDist>(new NoDist())};
  }
}




template<typename RNG>
void
fire(BrownionGSPN& et, TransitionType trans_id,
  BrownionState& s, smv::LocalMarking<Mark>& lm, RNG& rng)
{
  lm.template move<0,0>(0, 1, 1);
}



std::vector<std::tuple<place_t<BrownionGSPN>,size_t,int>>
neighbors_of_transition(BrownionGSPN& g, trans_t<BrownionGSPN> trans_id)
{
  std::vector<std::tuple<place_t<BrownionGSPN>,size_t,int>> place_ids;
  place_ids.push_back(std::make_tuple(trans_id.from, 0, -1));
  place_ids.push_back(std::make_tuple(trans_id.to, 0, -1));
  return place_ids;
}



template<typename F>
void neighbors_of_places(BrownionGSPN& g,
  const std::set<place_t<BrownionGSPN>>& place_id, const F& func)
{
  for (auto p : place_id)
  {
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

} // smv
} // afidd


#endif /* _BROWNION_MODEL_H_ */
