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
//   =========================================================================== 


#include <iostream>
#include <fstream>
#include <random>
#include <boost/cstdlib.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "smv.hpp"
#include "logging.hpp"


typedef std::size_t counter_type;
typedef std::string string_type;
typedef long coordinate_type;
typedef double real_type;

struct Walker
{
  // Intentionally empty
};

enum class Mobility : int { mobile = 1, immobile = -1 };

std::ostream& operator<<(std::ostream& os, const Mobility& m)
{
  return os << (m == Mobility::mobile ? "M" : "I");
};

struct PlaceKey
{
  PlaceKey() = default;
  PlaceKey(const PlaceKey&) = default;
  PlaceKey(coordinate_type _x, coordinate_type _y, Mobility _m)
    : x(_x), y(_y), mobility(_m) 
  {
    // Empty
  }

  coordinate_type x, y;
  Mobility mobility;

  friend inline 
  std::ostream& operator<<(std::ostream& os, const PlaceKey& pk) {
    return os << "{" << pk.x << ", " << pk.y << ", " << pk.mobility << "}";
  }

  friend inline 
  bool operator==(const PlaceKey& a, const PlaceKey& b) {
    return (a.x == b.x) && (a.y == b.y) && (a.mobility == b.mobility);
  }

  friend inline 
  bool operator<(const PlaceKey& a, const PlaceKey& b) {
    return afidd::smv::LazyLess(a.x, b.x, a.y, b.y, a.mobility, b.mobility);
  }
};

struct TransitionKey
{
  TransitionKey() = default;
  TransitionKey(const TransitionKey&) = default;
  TransitionKey(PlaceKey _from, PlaceKey _to) : 
    from(_from), to(_to)
  {
    // Empty
  }

  PlaceKey from, to;

  friend inline 
  std::ostream& operator<<(std::ostream& os, const TransitionKey& tk) {
    return os << "{" << tk.from << " -> " << tk.to << "}";
  }

  friend inline 
  bool operator==(const TransitionKey& a, const TransitionKey& b) {
    return (a.from == b.from) && (a.to == b.to);
  }

  friend inline 
  bool operator<(const TransitionKey& a, const TransitionKey& b) {
    return afidd::smv::LazyLess(a.from, b.from, a.to, b.to);
  }
};

using RandGen = std::mt19937;

using Dist    = afidd::smv::TransitionDistribution<RandGen>;
using ExpDist = afidd::smv::ExponentialDistribution<RandGen>;
using Weibull = afidd::smv::WeibullDistribution<RandGen>;

class BrownionGSPN
{
 public:
  typedef ::PlaceKey PlaceKey;
  typedef ::TransitionKey TransitionKey;
  // Could store the state parameters and distributions here
  // if we wanted.
};

struct UserState
{
  real_type weibull_shape, weibull_scale;

  void shape(const real_type& s) { weibull_shape=s; }
  void scale(const real_type& s) { weibull_scale=s; }

  const real_type& shape() const { return(weibull_shape); }
  const real_type& scale() const { return(weibull_scale); }
};


using TokenContainer = afidd::smv::Uncolored<Walker>;
using Mark  = afidd::smv::Marking<PlaceKey,
    TokenContainer>;
using Local = afidd::smv::LocalMarking<TokenContainer>;

std::pair<bool,std::unique_ptr<afidd::smv::TransitionDistribution<RandGen>>>
Enabled(const BrownionGSPN& et, TransitionKey trans_id,
        const UserState& s, const Local& lm, double te, double t0, RandGen& rng)
{
  if (lm.template Length<0>(0)>0) {
    // This is where we choose the distributions for the two
    // Brownion states.
    if (trans_id.from.mobility==Mobility::mobile) {
      return {true, std::unique_ptr<Weibull>(new Weibull(6.0, 3.0, te))};
    } else {
      return {true, std::unique_ptr<Dist>(new ExpDist(1.0, te))};
    }
  } else {
    return {false, std::unique_ptr<Dist>(nullptr)};
  }
}

void Fire(BrownionGSPN& et, TransitionKey trans_id, UserState& s, Local& lm,
      double t0, RandGen& rng)
{
  constexpr int token_layer = 0;
  const int from_edge = 0;
  const int to_edge = 1;
  const int cnt_to_move = 1;
  lm.template Move<token_layer,token_layer>(from_edge, to_edge, cnt_to_move);
}

std::vector<std::tuple<PlaceKey,int,int>>
NeighborsOfTransition(BrownionGSPN& g, TransitionKey trans_id)
{
  return {std::make_tuple(trans_id.from, 0, -1),
          std::make_tuple(trans_id.to, 0, 1)};
}

std::vector<std::tuple<PlaceKey,int,int>>
InputsOfTransition(BrownionGSPN& g, TransitionKey trans_id)
{
  return {std::make_tuple(trans_id.from, 0, -1),
          std::make_tuple(trans_id.to, 0, 1)};
}

template<typename F>
void NeighborsOfPlaces(BrownionGSPN& g,
    const std::set<PlaceKey>& place_id, const F& func)
{
  for (auto p : place_id) {
    // Transitions that start at this place.
    if (p.mobility == Mobility::mobile) {
      func(TransitionKey{p, {p.x, p.y-1, p.mobility}});
      func(TransitionKey{p, {p.x, p.y+1, p.mobility}});
      func(TransitionKey{p, {p.x-1, p.y, p.mobility}});
      func(TransitionKey{p, {p.x+1, p.y, p.mobility}});

      // Transitions that end at this place.
      func(TransitionKey{{p.x+1, p.y, p.mobility}, p});
      func(TransitionKey{{p.x-1, p.y, p.mobility}, p});
      func(TransitionKey{{p.x, p.y+1, p.mobility}, p});
      func(TransitionKey{{p.x, p.y-1, p.mobility}, p});

      // Transition to immobile
      func(TransitionKey{p, {p.x, p.y, Mobility::immobile}});
    } else {
      // Transition to mobile
      func(TransitionKey{p, {p.x, p.y, Mobility::mobile}});
    }
  }
}
