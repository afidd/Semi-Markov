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
/*!
 *  This example comes from a paper of partially-mixed metapopulations.
 *
 *  Hindes, Jason, Sarabjeet Singh, Christopher R. Myers, and David J.
 *  Schneider. "Epidemic fronts in complex networks with metapopulation
 *  structure." Physical Review E 88, no. 1 (2013): 012809.
 */


#include <memory>
#include <cmath>
#include <map>
#include <random>
#include <iostream>
#include "boost/program_options.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/binomial_distribution.hpp"
#include "boost/graph/erdos_renyi_generator.hpp"
#include "logging.hpp"
#include "smv.hpp"


using namespace afidd::smv;

using RNG=std::mt19937;


// The contact graph will take time to construct.
using ContactGraph=boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::undirectedS,
  boost::no_property,
  boost::no_property,
  boost::no_property,
  boost::listS
  >;


/*! Given a graph, find a list of nodes to whom to connect
 *  using a binomial distribution on the node degree.
 *  Return a list of nodes in random order.
 */
template<typename RNG>
std::vector<int64_t> AvailableToConnect(const ContactGraph& g,
    int64_t metapopulation_idx, int64_t nodes_per_metapopulation,
    double binomial_p, RNG& rng)
{
  int64_t offset=metapopulation_idx*nodes_per_metapopulation;

  std::vector<int64_t> connectable;
  for (int64_t node_idx=offset;
      node_idx<offset+nodes_per_metapopulation;
      ++node_idx)
  {
    auto degree=out_degree(node_idx, g);
    boost::random::binomial_distribution<> binomial(degree, binomial_p);
    auto connection_cnt=binomial(rng);
    for (int64_t connect_idx=0; connect_idx<connection_cnt; ++connect_idx) {
      connectable.push_back(node_idx);
    }
  }

  return connectable;
}




/*! Construct a sequence of erdos-renyi metapopulations connected
 *  with binomial samples on the edge degree of each node.
 */
template<typename RNG>
ContactGraph SerialMetapopulation(int64_t nodes_per_metapopulation,
    int64_t metapopulation_cnt, double edge_fraction, double binomial_p,
    RNG& rng)
{
  ContactGraph g(nodes_per_metapopulation*metapopulation_cnt);

  // Make multiple independent metapopulations, each from its own
  // erdos_renyi distribution.
  using ERGen=boost::erdos_renyi_iterator<RNG, ContactGraph>;
  for (int64_t make_meta_idx=0;
      make_meta_idx<metapopulation_cnt;
      ++make_meta_idx) {

    int64_t offset=make_meta_idx*nodes_per_metapopulation;

    ERGen erdos_renyi_iter(rng, nodes_per_metapopulation,
        edge_fraction);
    int64_t edge_cnt=0;
    for (; erdos_renyi_iter!=ERGen(); ++erdos_renyi_iter) {
      int64_t source_vertex, target_vertex;
      std::tie(source_vertex, target_vertex)=*erdos_renyi_iter;
      // remove self-loops here.
      if (source_vertex!=target_vertex) {
        add_edge(source_vertex+offset, target_vertex+offset, g);
        ++edge_cnt;
      }
    }
    BOOST_LOG_TRIVIAL(info)<<"metapop "<<make_meta_idx<<" has "
      <<edge_cnt<<" edges";
  }

  // Connect the distinct metapopulations.
  for (int64_t left_idx=0; left_idx<metapopulation_cnt-1; ++left_idx)
  {
    std::vector<int64_t> left_connectable=
      AvailableToConnect(g, left_idx, nodes_per_metapopulation,
      binomial_p, rng);
    std::vector<int64_t> right_connectable=
      AvailableToConnect(g, left_idx+1, nodes_per_metapopulation,
      binomial_p, rng);

    int64_t connection_cnt;
    if (left_connectable.size() > right_connectable.size()){
        std::shuffle(left_connectable.begin(), left_connectable.end(), rng);
        connection_cnt=right_connectable.size();
    } else {
        std::shuffle(right_connectable.begin(), right_connectable.end(), rng);
        connection_cnt=left_connectable.size();
    }

    for (int64_t cidx=0; cidx<connection_cnt; ++cidx) {
      add_edge(left_connectable[cidx], right_connectable[cidx], g);
    }
    BOOST_LOG_TRIVIAL(info)<<connection_cnt<<" edges between "<<left_idx
        <<" and "<<left_idx+1;
  }

  return g;
}




// Build the types for the model.
struct IndividualToken
{
  inline friend
  std::ostream& operator<<(std::ostream& os, const IndividualToken& it) {
    return os << 't';
  }
};



struct SIRPlace
{
  int disease;
  int64_t individual;
  int64_t metapop;

  SIRPlace()=default;
  SIRPlace(int d, int64_t i, int64_t m)
  : disease(d), individual(i), metapop(m)
  {}

  friend inline
  bool operator<(const SIRPlace& a, const SIRPlace& b) {
    return LazyLess(a.disease, b.disease, a.individual,
      b.individual, a.metapop, b.metapop);
  }


  friend inline
  bool operator==(const SIRPlace& a, const SIRPlace& b) {
    return (a.disease==b.disease)&& (a.individual==b.individual)
        && (a.metapop==b.metapop);
  }


  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRPlace& cp) {
    return os << '(' << cp.disease << ", " << cp.individual
        << ", "<<cp.metapop << ')';
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
  SIRPlace ind1;
  SIRPlace ind2;

  SIRTKey()=default;
  SIRTKey(SIRPlace c1, SIRPlace c2)
  : ind1(c1), ind2(c2)
  {}

  friend inline
  bool operator<(const SIRTKey& a, const SIRTKey& b) {
    return LazyLess(a.ind1, b.ind1, a.ind2, b.ind2);
  }

  friend inline
  bool operator==(const SIRTKey& a, const SIRTKey& b) {
    return (a.ind1==b.ind1) && (a.ind2==b.ind2);
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRTKey& cp) {
    return os << '(' << cp.ind1 << ", " << cp.ind2 << ')';
  }
};



using Local=LocalMarking<Uncolored<IndividualToken>>;
// Marking of the net.
using Mark=Marking<int64_t, Uncolored<IndividualToken>>;
// State of the continuous dynamical system.
class WithParams
{
  std::map<int,double> params;
};

using SIRGSPN=ExplicitTransitions<SIRPlace, SIRTKey, Local, RNG, WithParams>;
using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;
using Dist=TransitionDistribution<RNG>;
using ExpDist=ExponentialDistribution<RNG>;
using SIRTransition=ExplicitTransition<Local,RNG,WithParams>;


// Now make specific transitions.
class InfectNeighbor : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>> Enabled(
    const WithParams& s, const Local& lm, double te, double t0, RNG& rng) override {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(WithParams& s, Local& lm,
      double t0, RNG& rng) override {
    BOOST_LOG_TRIVIAL(debug) << "Fire infection " << lm;
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }

};





// Now make specific transitions.
class Recover : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>> Enabled(
      const WithParams& s, const Local& lm,
      double te, double t0, RNG& rng) override {
    if (lm.template InputTokensSufficient<0>()) {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(WithParams& s, Local& lm,
      double t0, RNG& rng) override {
    BOOST_LOG_TRIVIAL(debug) << "Fire recovery "<< lm;
    lm.template TransferByStochiometricCoefficient<0>(rng);
  }
};




SIRGSPN BuildSystem(int64_t nodes_per_metapopulation,
    int64_t metapopulation_cnt, double edge_fraction, double binomial_p,
    RNG& rng)
{
  auto contact_graph=SerialMetapopulation(nodes_per_metapopulation,
    metapopulation_cnt, edge_fraction, binomial_p, rng);

  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  enum { s, i, r };

  for (int64_t meta_idx=0; meta_idx<metapopulation_cnt; ++meta_idx) {
    for (int64_t ind_idx=0; ind_idx<nodes_per_metapopulation; ++ind_idx) {
      auto sp=SIRPlace{s, ind_idx, meta_idx};
      bg.AddPlace(sp);
      auto ip=SIRPlace{i, ind_idx, meta_idx};
      bg.AddPlace(ip);
      auto rp=SIRPlace{r, ind_idx, meta_idx};
      bg.AddPlace(rp);

      bg.AddTransition({ip, rp},
        {Edge{ip, -1}, Edge{rp, -1}},
        std::unique_ptr<SIRTransition>(new Recover()));
    }
  }

  auto index_to_place=[&](int64_t index)->SIRPlace {
    int64_t mpop=index / nodes_per_metapopulation;
    int64_t individual=index % nodes_per_metapopulation;
    return SIRPlace(s, individual, mpop);
  };

  std::set<std::array<int64_t,2>> seen;
  for (auto ce=edges(contact_graph); ce.first!=ce.second; ++ce.first) {
    int64_t left_index=static_cast<int64_t>(source(*ce.first, contact_graph));
    auto lefts=index_to_place(left_index);
    int64_t right_index=static_cast<int64_t>(target(*ce.first, contact_graph));
    auto rights=index_to_place(right_index);

    std::array<int64_t,2> new_transition{
      (std::min)(left_index, right_index), (std::max)(left_index, right_index)};
    if (seen.find(new_transition)==seen.end()) {
      SIRPlace lefti{i, lefts.individual, lefts.metapop};
      SIRPlace righti{i, rights.individual, rights.metapop};

      bg.AddTransition({lefti, rights},
        {Edge{lefti, -1}, Edge{rights, -1}, Edge{lefti, 1}, Edge{righti, 1}},
        std::unique_ptr<SIRTransition>(new InfectNeighbor()));

      bg.AddTransition({righti, lefts},
        {Edge{righti, -1}, Edge{lefts, -1}, Edge{righti, 1}, Edge{lefti, 1}},
        std::unique_ptr<SIRTransition>(new InfectNeighbor()));

      seen.insert(new_transition);
    }
  }

  return std::move(bg.Build());
}





template<typename GSPN>
class SIROutputFunction
{
public:
  typedef void result_type;

private:
  int64_t _pop_cnt;
  int64_t _ind_cnt;
  int64_t _threshold;
  std::vector<int64_t> _infected_per_pop;
  std::vector<double> _first_passage_time;
  std::vector<bool> _passed;

public:
  SIROutputFunction(int64_t population_cnt,
    int64_t individuals_per_metapopulation, int64_t threshold_for_passage)
  : _pop_cnt(population_cnt), _ind_cnt(individuals_per_metapopulation),
    _threshold(threshold_for_passage),
    _first_passage_time(population_cnt, 0), _passed(population_cnt, false),
    _infected_per_pop(population_cnt, 0)
  {}



  result_type operator()(const GSPN& gspn, const SIRState& state) {
    auto& modified=state.marking.Modified();
    for (auto place_idx : modified) {
      SIRPlace p=gspn.VertexPlace(place_idx);
      if (p.disease==1) {
        bool filled=(Length<0>(state.marking, place_idx)>0);
        if (filled) {
          _infected_per_pop[p.metapop]+=1;
          if (_infected_per_pop[p.metapop]==_threshold) {
            _first_passage_time[p.metapop]=state.CurrentTime();
            _passed[p.metapop]=true;
            BOOST_LOG_TRIVIAL(info)<<"Population "<<p.metapop
              <<" at time "<<state.CurrentTime();
          }
        }
      }
    }
  }
};





int main(int argc, char* argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Serial metapopulations");
  double poisson_constant=2.90156063;
  double binomial_to_neighbors=0.3;
  int64_t individuals_per_metapopulation=10000;
  int64_t metapopulation_cnt=20;
  int64_t rand_seed=1;
  std::string log_level;

  desc.add_options()
    ("help", "show help message")
    ("size,s",
      po::value<int64_t>(&individuals_per_metapopulation)->default_value(1000),
      "size of each partially-mixed population")
    ("seed,r",
      po::value<int64_t>(&rand_seed)->default_value(1),
      "seed for random number generator")
    ("count,c",
      po::value<int64_t>(&metapopulation_cnt)->default_value(20),
      "number of populations")
    ("erdos", po::value<double>(&poisson_constant)->default_value(2.90156063),
      "erdos-renyi constant")
    ("neighbor",
      po::value<double>(&binomial_to_neighbors)->default_value(0.3),
      "ratio of degree to neighbor given internal degree of node")
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

  afidd::LogInit(log_level);
  RNG rng(rand_seed);


  auto gspn=BuildSystem(individuals_per_metapopulation, metapopulation_cnt,
    poisson_constant/individuals_per_metapopulation, binomial_to_neighbors, rng);

  int64_t individual_cnt=metapopulation_cnt*individuals_per_metapopulation;
  BOOST_LOG_TRIVIAL(debug)<<"Initializing "<<individual_cnt<< " individuals "
      "as susceptible.";

  SIRState state;
  for (int64_t meta_idx=0; meta_idx<metapopulation_cnt; ++meta_idx) {
    for (int64_t ind_idx=0; ind_idx<individuals_per_metapopulation; ++ind_idx) {
      auto vert=gspn.PlaceVertex(SIRPlace{0, ind_idx, meta_idx});
      Add<0>(state.marking, vert, IndividualToken{});
    }
  }

  using Propagator=NonHomogeneousPoissonProcesses<decltype(gspn)::TransitionKey,
    RNG>;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RNG>;
  Propagator competing;
  Dynamics dynamics(gspn, {&competing});

  SIROutputFunction<SIRGSPN> output(
    metapopulation_cnt, individuals_per_metapopulation,
    std::lround(0.1*individuals_per_metapopulation));

  // The initial input string moves a token from susceptible to infected.
  auto first_case=static_cast<int64_t>(
      afidd::smv::uniform_index(rng, individuals_per_metapopulation));
  auto susceptible=gspn.PlaceVertex(SIRPlace{0, first_case, 0});
  auto infected=gspn.PlaceVertex(SIRPlace{1, first_case, 0});

  auto input_string=[&susceptible, &infected](SIRState& state)->void {
    Move<0,0>(state.marking, susceptible, infected, 1);
  };
  input_string(state);

  dynamics.Initialize(&state, &rng);

  int64_t transition_cnt=0;
  bool running=true;
  while (running) {
    output(gspn, state);
    running=dynamics(state);
    ++transition_cnt;
  }
  BOOST_LOG_TRIVIAL(info)<< transition_cnt << " transitions";
  return 0;
}

