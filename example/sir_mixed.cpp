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
/*! Fully mixed SIR with arbitrary transition rates.
 */

#include <tuple>
#include <map>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <memory>
#include <set>
#include <functional>
#include "stochnet.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/program_options.hpp"
#include "smv.hpp"

namespace smv=afidd::smv;
using namespace smv;
using RandGen=std::mt19937_64;


struct IndividualToken
{
  IndividualToken()=default;

  inline friend
  std::ostream& operator<<(std::ostream& os, const IndividualToken& it){
    return os << "T";
  }
};


struct SIRPlace
{
  int64_t disease;
  int64_t individual;

  SIRPlace()=default;
  SIRPlace(int64_t d, int64_t i) : disease(d), individual(i) {}

  friend inline
  bool operator<(const SIRPlace& a, const SIRPlace& b) {
    return LazyLess(a.disease, b.disease, a.individual,
      b.individual);
  }


  friend inline
  bool operator==(const SIRPlace& a, const SIRPlace& b) {
    return (a.disease==b.disease)&& (a.individual==b.individual);
  }


  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRPlace& cp) {
    return os << '(' << cp.disease << ", " << cp.individual<<')';
  }
};



struct SIRTKey
{
  int64_t ind1;
  int64_t ind2;
  int64_t kind;

  SIRTKey()=default;
  SIRTKey(int64_t c1, int64_t c2, int64_t k) : ind1(c1), ind2(c2), kind(k) {}

  friend inline
  bool operator<(const SIRTKey& a, const SIRTKey& b) {
    return LazyLess(a.ind1, b.ind1, a.ind2, b.ind2,
      a.kind, b.kind);
  }

  friend inline
  bool operator==(const SIRTKey& a, const SIRTKey& b) {
    return (a.ind1==b.ind1) && (a.ind2==b.ind2) && (a.kind==b.kind);
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRTKey& cp) {
    return os << '(' << cp.ind1 << ", " << cp.ind2 << ", "
     << cp.kind << ')';
  }
};



// This is as much of the marking as the transition will see.
using Local=LocalMarking<Uncolored<IndividualToken>>;
// Extra state to add to the system state. Will be passed to transitions.
struct WithParams {
  // Put our parameters here.
  std::map<int,double> params;
};


// The transition needs to know the local marking and any extra state.
using SIRTransition=ExplicitTransition<Local,RandGen,WithParams>;

using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;



// Now make specific transitions.
class InfectNeighbor : public SIRTransition
{
public:
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(s.params.at(0), te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};





// Now make specific transitions.
class Recover : public SIRTransition
{
public:
  virtual
  std::pair<bool, std::unique_ptr<Dist>> Enabled(const UserState& s,
      const Local& lm, double te, double t0) const override {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(s.params.at(1), te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire recovery "<< lm;
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};


// The GSPN itself.
using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;

/*! SIR infection on an all-to-all graph of uncolored tokens.
 */
SIRGSPN
BuildSystem(int64_t individual_cnt)
{
  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  enum { s, i, r };

  for (int64_t ind_idx=0; ind_idx<individual_cnt; ind_idx++) {
    for (int64_t place : std::vector<int>{s, i, r}) {
      bg.AddPlace({place, ind_idx}, 0);
    }
  }

  for (int64_t left_idx=0; left_idx<individual_cnt-1; left_idx++) {
    bg.AddTransition({left_idx, left_idx, 0},
      {Edge{{i, left_idx}, -1}, Edge{{r, left_idx}, 1}},
      std::unique_ptr<SIRTransition>(new Recover())
      );

    for (int64_t right_idx=left_idx+1; right_idx<individual_cnt; right_idx++) {
      SIRPlace left{i, left_idx};
      SIRPlace rights{s, right_idx};
      SIRPlace righti{i, right_idx};

      bg.AddTransition({left_idx, right_idx, 0},
        {Edge{left, -1}, Edge{rights, -1}, Edge{left, 1}, Edge{righti, 1}},
        std::unique_ptr<SIRTransition>(new InfectNeighbor()));

      SIRPlace lefts{s, left_idx};
      SIRPlace lefti{i, left_idx};
      SIRPlace right{i, right_idx};

      bg.AddTransition({right_idx, left_idx, 0},
        {Edge{right, -1}, Edge{lefts, -1}, Edge{right, 1}, Edge{lefti, 1}},
        std::unique_ptr<SIRTransition>(new InfectNeighbor()));
    }
  }

  // std::move the transitions because they contain unique_ptr.
  return std::move(bg.Build());
}




/*! Write a file showing place ids and the internal indices for debugging.
 */
template<typename GSPN>
void WriteIds(const GSPN& gspn, const std::string& fname,
    int64_t individual_cnt)
{
  std::ofstream out{fname};

  for (int64_t individual=0; individual<individual_cnt; ++individual) {
    for (int64_t disease=0; disease<3; disease++) {
      SIRPlace p{disease, individual};
      auto vert=gspn.PlaceVertex(p);
      out << vert << '\t' << p << std::endl;
    }
  }
}


template<typename SIRState, typename SIRGSPN>
struct SIROutput
{
  int64_t step_cnt{0};
  const SIRGSPN& gspn_;

  SIROutput(const SIRGSPN& gspn) : gspn_(gspn) {}

  void operator()(const SIRState& state) {
    ++step_cnt;
    BOOST_LOG_TRIVIAL(debug)<<"trans "<<state.last_transition;
    auto transition_key=gspn_.VertexTransition(state.last_transition);
    int64_t id0, id1;
    typename SIRGSPN::UserPlaceKey key0, key1;

    std::tie(id0, key0)=gspn_.PlaceOfTransition(state.last_transition, 0);
    std::tie(id1, key1)=gspn_.PlaceOfTransition(state.last_transition, 1);
    BOOST_LOG_TRIVIAL(debug)<<"place0 "<<key0<<" place1 "<<key1;
    BOOST_LOG_TRIVIAL(debug) << "trans " << transition_key
      << " time " << state.CurrentTime() << " step " << step_cnt;
    BOOST_LOG_TRIVIAL(trace) << state.marking;
  }

  void final(const SIRState& state) {
    BOOST_LOG_TRIVIAL(info) << "Took "<< step_cnt << " transitions.";
  }
};



int main(int argc, char *argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR");
  int64_t individual_cnt=10;
  size_t rand_seed=1;
  double beta=1.0;
  double gamma=1.0;
  std::string log_level;
  std::string translation_file;

  desc.add_options()
    ("help", "show help message")
    ("size,s",
      po::value<int64_t>(&individual_cnt)->default_value(10),
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
    ("translate",
      po::value<std::string>(&translation_file)->default_value(""),
      "write file relating place ids to internal ids")
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

  afidd::LogInit(log_level);
  RandGen rng(rand_seed);

  auto gspn=BuildSystem(individual_cnt);

  if (translation_file.size()>0)
  {
    WriteIds(gspn, translation_file, individual_cnt);
  }


  // Marking of the net.
  static_assert(std::is_same<int64_t,SIRGSPN::PlaceKey>::value,
    "The GSPN's internal place type is int64_t.");
  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<IndividualToken>>;
  using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

  SIRState state;
  state.user.params[0]=beta;
  state.user.params[1]=gamma;

  for (int64_t individual=0; individual<individual_cnt; ++individual) {
    auto susceptible=gspn.PlaceVertex({0, individual});
    Add<0>(state.marking, susceptible, IndividualToken{});
  }

  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  BOOST_LOG_TRIVIAL(debug) << state.marking;

  // The initial input string moves a token from susceptible to infected.
  auto first_case=static_cast<int64_t>(
      smv::uniform_index(rng, individual_cnt));
  BOOST_LOG_TRIVIAL(trace)<<"First case is "<<first_case;
  int64_t first_s=gspn.PlaceVertex({0, first_case});
  int64_t first_i=gspn.PlaceVertex({1, first_case});
  auto input_string=[&first_s, &first_i](SIRState& state)->void {
    Move<0,0>(state.marking, first_s, first_i, 1);
  };

  input_string(state);

  SIROutput<SIRState,SIRGSPN> output_function(gspn);

  dynamics.Initialize(&state, &rng);

  bool running=true;
  auto nothing=[](SIRState&)->void {};
  while (running) {
    running=dynamics(state);
    if (running) output_function(state);
  }
  output_function.final(state);
  return 0;
}

