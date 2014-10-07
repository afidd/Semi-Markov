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
#include <mutex>
#include "stochnet.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "ensemble.hpp"

namespace smv=afidd::smv;
using namespace smv;
using RandGen=std::mt19937_64;


enum class SIRParam { Beta, Gamma };
/*! This class takes a parameter from the command line or a file,
 *  with enough information to describe and store it.
 */
struct Parameter {
  SIRParam kind;
  std::string name;
  double value;
  std::string description;
  Parameter(SIRParam k, std::string n, double v, std::string desc)
  : kind(k), name(n), value(v), description(desc) {}
  Parameter()=default;
};

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
  std::map<SIRParam,double> params;
  int64_t token_cnt;
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
    double te, double t0, RandGen& rng) override {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(s.params.at(SIRParam::Beta), te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    //lm.template TransferByStochiometricCoefficient<0>(rng);
    lm.template Move<0,0>(1, 3, 1);
  }
};





// Now make specific transitions.
class Recover : public SIRTransition
{
public:
  virtual
  std::pair<bool, std::unique_ptr<Dist>> Enabled(const UserState& s,
      const Local& lm, double te, double t0, RandGen& rng) override {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(s.params.at(SIRParam::Gamma), te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    BOOST_LOG_TRIVIAL(trace) << "Fire recovery "<< lm;
    //lm.template TransferByStochiometricCoefficient<0>(rng);
    lm.template Move<0,0>(0, 1, 1);
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

/*! An observer sees the whole state, but this is the part we extract.
 *  This struct is what will be pulled from the state at each
 *  time step.
 */
struct TrajectoryEntry {
  int64_t s;
  int64_t i;
  int64_t r;
  double t;
  TrajectoryEntry(int64_t s, int64_t i, int64_t r, double t)
  : s(s), i(i), r(r), t(t) {}
  TrajectoryEntry()=default;
};

/*! Abstract base for observers of the trajectory,
 *  because maybe we just want to save the final recovered count,
 *  but maybe we want the area under the curve.
 */
class TrajectoryObserver
{
public:
  virtual void Step(TrajectoryEntry sirt)=0;
  virtual const std::vector<TrajectoryEntry>& Trajectory() const =0;
};


/*! Save the whole trajectory.
 */
class TrajectorySave : public TrajectoryObserver
{
  std::vector<TrajectoryEntry> trajectory_;
  int64_t cnt_;
 public:
  TrajectorySave() : trajectory_{10000}, cnt_{0} {}
  virtual void Step(TrajectoryEntry sirt) override {
    trajectory_[cnt_]=sirt;
    ++cnt_;
    if (cnt_==trajectory_.size()) {
      // write everything
      cnt_=0;
    }
  }
  virtual const std::vector<TrajectoryEntry>& Trajectory() const {
    return trajectory_; }
  int64_t final_removed() const { return trajectory_[cnt_-1].r; }
  friend
  std::ostream& operator<<(std::ostream& os, const TrajectorySave& ts) {
    for (int64_t i=0; i<ts.cnt_; ++i) {
      os << ts.trajectory_[i].s << '\t' << ts.trajectory_[i].i
          << '\t' << ts.trajectory_[i].r << '\t'
          << ts.trajectory_[i].t << std::endl;
    }
    return os;
  }
};


template<typename SIRState, typename SIRGSPN>
struct SIROutput
{
  int64_t step_cnt{0};
  std::vector<int64_t> sir_;
  const SIRGSPN& gspn_;
  std::shared_ptr<TrajectoryObserver> observer_;

  SIROutput(const SIRGSPN& gspn, std::vector<int64_t> sir,
    std::shared_ptr<TrajectoryObserver> observer)
    : gspn_(gspn), sir_{sir}, observer_{observer} {}

  void operator()(const SIRState& state) {
    BOOST_LOG_TRIVIAL(debug)<<"trans "<<state.last_transition;
    auto transition_key=gspn_.VertexTransition(state.last_transition);
    if (transition_key.ind1==transition_key.ind2) {
      sir_[1]-=1;
      sir_[2]+=1;
    } else {
      sir_[0]-=1;
      sir_[1]+=1;
    }
    int64_t id0, id1;
    typename SIRGSPN::UserPlaceKey key0, key1;

    std::tie(id0, key0)=gspn_.PlaceOfTransition(state.last_transition, 0);
    std::tie(id1, key1)=gspn_.PlaceOfTransition(state.last_transition, 1);
    BOOST_LOG_TRIVIAL(debug)<<"place0 "<<key0<<" place1 "<<key1;
    BOOST_LOG_TRIVIAL(debug) << "trans " << transition_key
      << " time " << state.CurrentTime() << " step " << step_cnt;
    BOOST_LOG_TRIVIAL(trace) << state.marking;

    //assert(this->check_sir(state));
    observer_->Step({sir_[0], sir_[1], sir_[2], state.CurrentTime()});
    ++step_cnt;
  }

  bool check_sir(const SIRState& state) {
    int64_t sir_cnt=sir_[0]+sir_[1]+sir_[2];
    int64_t token_id=0;
    std::vector<int64_t> sir(3);
    for (int64_t sir_idx=0; sir_idx<3; ++sir_idx) {
      sir[sir_idx]=0;
      for (int64_t sus_idx=0; sus_idx<sir_cnt; ++sus_idx) {
        sir[sir_idx]+=Length<0>(state.marking,
            gspn_.PlaceVertex(SIRPlace{sir_idx, sus_idx}));
      }
    }
    for (int cidx=0; cidx<sir.size(); ++cidx) {
      if (sir[cidx]!=sir_[cidx]) {
        BOOST_LOG_TRIVIAL(error)<< "Expected: "
            <<sir_[0]<<'\t'<<sir_[1]<<'\t'<<sir_[2];
        BOOST_LOG_TRIVIAL(error)<< "Found: "
            <<sir[0]<<'\t'<<sir[1]<<'\t'<<sir[2];
        return false;
      }
    }
    return true;
  }

  void final(const SIRState& state) {
  }
};


int64_t SIR_run(const std::vector<int64_t>& sir_cnt, SIRGSPN& gspn,
    const std::vector<Parameter>& parameters,
    std::shared_ptr<TrajectoryObserver> observer,
    RandGen& rng) {
  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<IndividualToken>>;
  using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

  SIRState state;
  for (auto& cp : parameters) {
    state.user.params[cp.kind]=cp.value;
  }

  int64_t token_id=0;
  for (int64_t sir_idx=0; sir_idx<3; ++sir_idx) {
    for (int64_t sus_idx=0; sus_idx<sir_cnt[sir_idx]; ++sus_idx) {
      Add<0>(state.marking, gspn.PlaceVertex(SIRPlace{sir_idx, sus_idx}),
        IndividualToken{});
      ++token_id;
    }
  }
  state.user.token_cnt=token_id;

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  PropagateCompetingProcesses<int64_t,RandGen> simple;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing, &simple});

  BOOST_LOG_TRIVIAL(debug) << state.marking;

  SIROutput<SIRState,SIRGSPN> output_function(gspn, sir_cnt, observer);

  dynamics.Initialize(&state, &rng);

  bool running=true;
  auto nothing=[](SIRState&)->void {};
  double last_time=state.CurrentTime();
  while (running) {
    SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"SIR_run() time "<<state.CurrentTime());
    running=dynamics(state);
    if (running) {
      double new_time=state.CurrentTime();
      if (new_time-last_time<0) {
        BOOST_LOG_TRIVIAL(warning) << "last time "<<last_time <<" "
          << " new_time "<<new_time;
      }
      last_time=new_time;
      output_function(state);
    }
  }
  if (running) {
    BOOST_LOG_TRIVIAL(info)<<"Reached end time "<<state.CurrentTime();
  } else {
    BOOST_LOG_TRIVIAL(info)<<"No transitions left to fire at time "<<last_time;
  }
  output_function.final(state);
  return 0;
}




int main(int argc, char *argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR");
  int64_t individual_cnt=10;
  int64_t thread_cnt=1;
  int64_t run_cnt=1;
  size_t rand_seed=1;
  std::string log_level;
  std::string translation_file;

  std::vector<Parameter> parameters;
  parameters.emplace_back(Parameter{SIRParam::Beta, "beta", 1,
    "main infection rate"});
  parameters.emplace_back(Parameter{SIRParam::Gamma, "gamma", 1,
    "recovery rate"});

  desc.add_options()
    ("help", "show help message")
    ("size,s",
      po::value<int64_t>(&individual_cnt)->default_value(10),
      "size of the population")
    ("threadcnt,j",
      po::value<int64_t>(&thread_cnt)->default_value(1),
      "number of threads to use")
    ("run",
      po::value<int64_t>(&run_cnt)->default_value(1),
      "number of runs of the model")
    ("seed,r",
      po::value<size_t>(&rand_seed)->default_value(1),
      "seed for random number generator")
    ("loglevel", po::value<std::string>(&log_level)->default_value("info"),
      "Set the logging level to trace, debug, info, warning, error, or fatal.")
    ("translate",
      po::value<std::string>(&translation_file)->default_value(""),
      "write file relating place ids to internal ids")
    ;

  for (auto& p : parameters) {
    desc.add_options()(p.name.c_str(),
      po::value<double>(&p.value)->default_value(p.value),
      p.description.c_str());
  }

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

  auto gspn=BuildSystem(individual_cnt);

  if (translation_file.size()>0)
  {
    WriteIds(gspn, translation_file, individual_cnt);
  }
  // Marking of the net.
  static_assert(std::is_same<int64_t,SIRGSPN::PlaceKey>::value,
    "The GSPN's internal place type is int64_t.");

  std::vector<int64_t> sir_init{individual_cnt-1, 1, 0};
  BOOST_LOG_TRIVIAL(info)<<"Starting with sir="<<sir_init[0]<<" "<<sir_init[1]
    <<" "<<sir_init[2];

  std::mutex printing;
  auto runnable=[=, &gspn, &printing]
      (RandGen& rng, size_t single_seed, size_t idx)->void {
    std::shared_ptr<TrajectorySave> observer=std::make_shared<TrajectorySave>();
    SIR_run(sir_init, gspn, parameters, observer, rng);
    std::unique_lock<std::mutex> haveprint(printing);
    std::cout << observer->final_removed() << std::endl;
    //std::cout << *observer;
    //file.SaveTrajectory(parameters, single_seed, idx, observer->Trajectory());
  };

  Ensemble<decltype(runnable),RandGen> ensemble(runnable, thread_cnt,
      run_cnt, rand_seed);
  ensemble.Run();
  BOOST_LOG_TRIVIAL(debug)<<"Finished running ensemble.";

  //file.Close();


  return 0;
}

