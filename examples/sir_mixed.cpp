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
#include "boost/program_options.hpp"
#include "gspn.h"
#include "petri_graph.h"
#include "marking.h"
#include "distributions.h"
#include "continuous_state.h"
#include "explicit_transitions.h"
#include "partial_core_matrix.h"
#include "continuous_dynamics.h"
#include "build_graph.h"
#include "smv_algorithm.h"
#include "logging.h"

namespace smv=afidd::smv;
using namespace smv;
using RandGen=std::mt19937_64;


struct IndividualToken
{
  IndividualToken()=default;

  inline friend
  std::ostream& operator<<(std::ostream& os, const IndividualToken& it)
  {
    return os << "T";
  }
};


struct SIRPlace
{
  size_t disease;
  size_t individual;

  SIRPlace()=default;
  SIRPlace(size_t d, size_t i)
  : disease(d), individual(i)
  {}

  friend inline
  bool operator<(const SIRPlace& a, const SIRPlace& b)
  {
    return lazy_less(a.disease, b.disease, a.individual,
      b.individual);
  }


  friend inline
  bool operator==(const SIRPlace& a, const SIRPlace& b)
  {
    return (a.disease==b.disease)&& (a.individual==b.individual);
  }


  friend inline
  std::ostream&
  operator<<(std::ostream& os, const SIRPlace& cp)
  {
    return os << '(' << cp.disease << ", " << cp.individual<<')';
  }
};



/*! This identifies a cow transition.
 *  We are being luxurious. If it's a cow-to-cow infection,
 *  both cows identify the transition. Subgroup-to-subgroup
 *  can also be recorded. The kind is then an identifier for
 *  a particular infection or movement.
 */
struct SIRTKey
{
  size_t ind1;
  size_t ind2;
  size_t kind;

  SIRTKey()=default;
  SIRTKey(size_t c1, size_t c2, size_t k)
  : ind1(c1), ind2(c2), kind(k)
  {}

  friend inline
  bool operator<(const SIRTKey& a, const SIRTKey& b)
  {
    return lazy_less(a.ind1, b.ind1, a.ind2, b.ind2,
      a.kind, b.kind);
  }

  friend inline
  bool operator==(const SIRTKey& a, const SIRTKey& b)
  {
    return (a.ind1==b.ind1) && (a.ind2==b.ind2) && (a.kind==b.kind);
  }

  friend inline
  std::ostream&
  operator<<(std::ostream& os, const SIRTKey& cp)
  {
    return os << '(' << cp.ind1 << ", " << cp.ind2 << ", "
     << cp.kind << ')';
  }
};


// Marking of the net.
using Mark=Marking<size_t, Uncolored<IndividualToken>>;
// State of the continuous dynamical system.
template<typename Mark, typename TimeStrategy=KahanTime>
class WithParamsState : public TimeStrategy
{
public:
  typedef Mark Marking;

  Marking marking;
  std::map<int,double> params;
};

using SIRState=WithParamsState<Mark>;


// This class holds the transitions.
using SIRGSPN=
    ExplicitTransitions<SIRState, SIRPlace, SIRTKey, RandGen>;


class SIRTransition
: public SIRGSPN::Transition
{
public:
  SIRTransition() {}
  virtual ~SIRTransition() {}
};



using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;
using NoDist=NoDistribution<RandGen>;




// Now make specific transitions.
class InfectNeighbor : public SIRTransition
{

  virtual std::pair<bool, std::unique_ptr<Dist>>
  enabled(const SIRState& s, const LocalMarking<Mark>& lm, double te) const override
  {
    if (lm.template input_tokens_sufficient<0>())
    {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(s.params.at(0), te))};
    }
    else
    {
      return {false, std::unique_ptr<NoDist>(new NoDist())};
    }
  }

  virtual void fire(SIRState& s, LocalMarking<Mark>& lm,
      RandGen& rng) const override
  {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template transfer_by_stochiometric_coefficient<0>(rng);
  }

};





// Now make specific transitions.
class Recover : public SIRTransition
{

  virtual std::pair<bool, std::unique_ptr<Dist>>
  enabled(const SIRState& s, const LocalMarking<Mark>& lm, double te) const override
  {
    if (lm.template input_tokens_sufficient<0>())
    {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(s.params.at(1), te))};
    }
    else
    {
      return {false, std::unique_ptr<NoDist>(new NoDist())};
    }
  }

  virtual void fire(SIRState& s, LocalMarking<Mark>& lm,
      RandGen& rng) const override
  {
    BOOST_LOG_TRIVIAL(trace) << "Fire recovery "<< lm;
    lm.template transfer_by_stochiometric_coefficient<0>(rng);
  }

};



/*! SIR infection on an all-to-all graph of uncolored tokens.
 */
SIRGSPN
build_system(size_t individual_cnt)
{
  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  enum { s, i, r };

  for (size_t ind_idx=0; ind_idx<individual_cnt; ind_idx++)
  {
    for (size_t place : std::vector<int>{s, i, r})
    {
      bg.add_place({place, ind_idx}, 0);
    }
  }

  for (size_t left_idx=0; left_idx<individual_cnt-1; left_idx++)
  {
    bg.add_transition({left_idx, left_idx, 0},
      {Edge{{i, left_idx}, -1}, Edge{{r, left_idx}, 1}},
      std::unique_ptr<SIRTransition>(new Recover())
      );

    for (size_t right_idx=left_idx+1; right_idx<individual_cnt; right_idx++)
    {
      SIRPlace left{i, left_idx};
      SIRPlace rights{s, right_idx};
      SIRPlace righti{i, right_idx};

      bg.add_transition({left_idx, right_idx, 0},
        {Edge{left, -1}, Edge{rights, -1}, Edge{left, 1}, Edge{righti, 1}},
        std::unique_ptr<SIRTransition>(new InfectNeighbor()));

      SIRPlace lefts{s, left_idx};
      SIRPlace lefti{i, left_idx};
      SIRPlace right{i, right_idx};

      bg.add_transition({right_idx, left_idx, 0},
        {Edge{right, -1}, Edge{lefts, -1}, Edge{right, 1}, Edge{lefti, 1}},
        std::unique_ptr<SIRTransition>(new InfectNeighbor()));
    }
  }

  // std::move the transitions because they contain unique_ptr.
  return std::move(bg.build());
}





int main(int argc, char *argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR");
  size_t individual_cnt=10;
  size_t rand_seed=1;
  double beta=1.0;
  double gamma=1.0;
  std::string log_level;

  desc.add_options()
    ("help", "show help message")
    ("size,s",
      po::value<size_t>(&individual_cnt)->default_value(10),
      "size of the population")
    ("seed,r",
      po::value<size_t>(&rand_seed)->default_value(1),
      "seed for random number generator")
    ("beta",
      po::value<double>(&beta)->default_value(1.0),
      "parameter for infection of neighbor")
    ("gamma",
      po::value<double>(&gamma)->default_value(1.0),
      "parameter for recovery")
    ("loglevel", po::value<std::string>(&log_level)->default_value("info"),
      "Set the logging level to trace, debug, info, warning, error, or fatal.")
    ;

  po::variables_map vm;
  auto parsed_options=po::parse_command_line(argc, argv, desc);
  po::store(parsed_options, vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    return 0;
  }

  afidd::log_init(log_level);
  RandGen rng(rand_seed);

  auto gspn=build_system(individual_cnt);

  SIRState state;
  state.params[0]=beta;
  state.params[1]=gamma;

  for (size_t individual=0; individual<individual_cnt; ++individual)
  {
    auto susceptible=gspn.place_vertex({0, individual});
    add<0>(state.marking, susceptible, IndividualToken{});
  }

  using Markov=PartialCoreMatrix<SIRGSPN, SIRState, RandGen>;
  Markov system(gspn, state);

  BOOST_LOG_TRIVIAL(debug) << state.marking;

  // The initial input string moves a token from susceptible to infected.
  auto first_case=smv::uniform_index(rng, individual_cnt);
  size_t first_s=gspn.place_vertex({0, first_case});
  size_t first_i=gspn.place_vertex({1, first_case});
  auto input_string=[&first_s, &first_i](SIRState& state)->void {
    move<0,0>(state.marking, first_s, first_i, 1);
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

