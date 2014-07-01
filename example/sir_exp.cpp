
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
#include "boost/math/constants/constants.hpp"
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
  int disease;
  SIRPlace()=default;
  SIRPlace(int d) : disease(d) {}
  friend inline
  bool operator<(const SIRPlace& a, const SIRPlace& b) {
    return LazyLess(a.disease, b.disease);
  }


  friend inline
  bool operator==(const SIRPlace& a, const SIRPlace& b) {
    return a.disease==b.disease;
  }


  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRPlace& cp) {
    return os << '(' << cp.disease << ')';
  }
};


struct SIRTKey
{
  int kind;
  SIRTKey()=default;
  SIRTKey(int k) : kind(k) {}

  friend inline
  bool operator<(const SIRTKey& a, const SIRTKey& b) {
    return LazyLess(a.kind, b.kind);
  }

  friend inline
  bool operator==(const SIRTKey& a, const SIRTKey& b) {
    return a.kind==b.kind;
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRTKey& cp) {
    return os << '(' << cp.kind << ')';
  }
};


enum class SIRParam { Beta0, Beta1, Gamma, Birth, Mu };

// This is as much of the marking as the transition will see.
using Local=LocalMarking<Uncolored<IndividualToken>>;
// Extra state to add to the system state. Will be passed to transitions.
struct WithParams {
  // Put our parameters here.
  std::map<SIRParam,double> params;
};


// The transition needs to know the local marking and any extra state.
using SIRTransition=ExplicitTransition<Local,RandGen,WithParams>;

using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;



// Now make specific transitions.
class Infect : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    auto S=lm.template Length<0>(0);
    auto I=lm.template Length<0>(1);
    auto R=lm.template Length<0>(2);
    if (S>0 && I>0) {
      double rate=S*s.params.at(SIRParam::Beta0)*
        (1+s.params.at(SIRParam::Beta0)*
          std::cos(2*boost::math::constants::pi<double>()*t0))/
          (S+I+R);
      return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template Move<0,0>(0, 3, 1);
  }
};



// Now make specific transitions.
class Recover : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    auto I=lm.template Length<0>(0);
    if (I>0) {
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(I*s.params.at(SIRParam::Gamma), te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template Move<0, 0>(0, 1, 1);
  }
};



// Now make specific transitions.
class Birth : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    return {true, std::unique_ptr<ExpDist>(
      new ExpDist(s.params.at(SIRParam::Birth), te))};
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template Add<0>(1, IndividualToken{});
  }
};



// Now make specific transitions.
class Death : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    auto SIR=lm.template Length<0>(0);
    if (SIR>0) {
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(SIR*s.params.at(SIRParam::Mu), te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template Remove<0>(0, 1, rng);
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

  for (int place : std::vector<int>{s, i, r}) {
    bg.AddPlace({place}, 0);
  }

  enum { infect, recover, birth, deaths, deathi, deathr };

  bg.AddTransition({infect},
    {Edge{{s}, -1}, Edge{{i}, -1}, Edge{{r}, -1}, Edge{{i}, 2}, Edge{{r}, 1}},
    std::unique_ptr<SIRTransition>(new Infect())
    );

  bg.AddTransition({recover},
    {Edge{{i}, -1}, Edge{{r}, 1}},
    std::unique_ptr<SIRTransition>(new Recover())
    );

  bg.AddTransition({birth},
    {Edge{{s}, -1}, Edge{{s}, 2}},
    std::unique_ptr<SIRTransition>(new Birth())
    );

  bg.AddTransition({deaths},
    {Edge{{s}, -1}, Edge{{s}, 0}},
    std::unique_ptr<SIRTransition>(new Death())
    );

  bg.AddTransition({deathi},
    {Edge{{i}, -1}, Edge{{i}, 0}},
    std::unique_ptr<SIRTransition>(new Death())
    );

  bg.AddTransition({deathr},
    {Edge{{r}, -1}, Edge{{r}, 0}},
    std::unique_ptr<SIRTransition>(new Death())
    );

  // std::move the transitions because they contain unique_ptr.
  return std::move(bg.Build());
}



template<typename SIRState>
struct SIROutput
{
  std::vector<int64_t> places_;
  using StateArray=std::array<int64_t,3>;
  std::vector<std::tuple<StateArray,double>> trajectory_;
  double max_time_;
  int64_t max_count_;

  SIROutput(double max_time, int64_t max_count,
    const std::vector<int64_t>& sir_places)
  : places_{sir_places},
    max_time_(max_time), max_count_(max_count)
  {};

  int64_t step_cnt{0};

  void operator()(const SIRState& state) {
    ++step_cnt;
    auto S=Length<0>(state.marking, places_[0]);
    auto I=Length<0>(state.marking, places_[1]);
    auto R=Length<0>(state.marking, places_[2]);
    //std::cout << "(" << S << "," << I << "," << R << ") "
    //  << state.CurrentTime() << std::endl;
    trajectory_.emplace_back(std::make_tuple(
      StateArray{S, I, R}, state.CurrentTime()));

  }

  void final(const SIRState& state) {
    BOOST_LOG_TRIVIAL(info) << "Took "<< step_cnt << " transitions.";
  }
};



int main(int argc, char *argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR");
  int64_t individual_cnt=100000;
  int64_t infected_start=std::floor(individual_cnt*0.001);
  int64_t recovered_start=std::floor(individual_cnt*0.9);
  size_t rand_seed=1;
  // Time is in years.
  double beta0=400; // Rate of infection is over one per day.
  double beta1=0.6;
  double gamma=365/14.0; // Rate of recovery is a rate per year.
  double deathrate=1/70.0;
  double birthrate=deathrate;
  double end_time=30.0;
  std::string log_level;
  std::string translation_file;

  desc.add_options()
    ("help", "show help message")
    ("size,s",
      po::value<int64_t>(&individual_cnt)->default_value(individual_cnt),
      "size of the population")
    ("seed,r",
      po::value<size_t>(&rand_seed)->default_value(rand_seed),
      "seed for random number generator")
    ("beta0",
      po::value<double>(&beta0)->default_value(beta0),
      "parameter beta0 for infection of neighbor")
    ("beta1",
      po::value<double>(&beta1)->default_value(beta1),
      "parameter beta1 for infection of neighbor")
    ("gamma",
      po::value<double>(&gamma)->default_value(gamma),
      "parameter for recovery")
    ("death",
      po::value<double>(&deathrate)->default_value(deathrate),
      "parameter for death")
    ("birth",
      po::value<double>(&birthrate)->default_value(birthrate),
      "parameter for birth")
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
  }


  // Marking of the net.
  static_assert(std::is_same<int64_t,SIRGSPN::PlaceKey>::value,
    "The GSPN's internal place type is int64_t.");
  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<IndividualToken>>;
  using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

  SIRState state;
  state.user.params[SIRParam::Beta0]=beta0;
  state.user.params[SIRParam::Beta1]=beta1;
  state.user.params[SIRParam::Gamma]=gamma;
  // Birthrate is not frequency-dependent. It scales differently
  // which creates a fixed point in the phase plane.
  state.user.params[SIRParam::Birth]=birthrate*individual_cnt;
  state.user.params[SIRParam::Mu]=deathrate;

  auto susceptible_place=gspn.PlaceVertex({0});
  auto susc_start=individual_cnt-infected_start-recovered_start;
  for (int64_t sus_idx=0; sus_idx<susc_start; ++sus_idx) {
    Add<0>(state.marking, susceptible_place, IndividualToken{});
  }
  auto infected_place=gspn.PlaceVertex({1});
  for (int64_t inf_idx=0; inf_idx<infected_start; ++inf_idx) {
    Add<0>(state.marking, infected_place, IndividualToken{});
  }
  auto recovered_place=gspn.PlaceVertex({2});
  for (int64_t rec_idx=0; rec_idx<recovered_start; ++rec_idx) {
    Add<0>(state.marking, recovered_place, IndividualToken{});
  }

  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  BOOST_LOG_TRIVIAL(debug) << state.marking;

  SIROutput<SIRState> output_function(end_time, individual_cnt*2,
    {susceptible_place, infected_place,
      recovered_place});

  dynamics.Initialize(&state, &rng);

  bool running=true;
  auto nothing=[](SIRState&)->void {};
  while (running && state.CurrentTime()<end_time) {
    running=dynamics(state);
    output_function(state);
  }
  output_function.final(state);
  return 0;
}

