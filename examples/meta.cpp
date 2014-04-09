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
#include "petri_graph.hpp"
#include "gspn.hpp"
#include "continuous_state.hpp"
#include "marking.hpp"
#include "distributions.hpp"
#include "explicit_transitions.hpp"
#include "partial_core_matrix.hpp"
#include "continuous_dynamics.hpp"
#include "smv_algorithm.hpp"
#include "local_marking.hpp"


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
std::vector<size_t> available_to_connect(const ContactGraph& g,
    size_t metapopulation_idx, size_t nodes_per_metapopulation,
    double binomial_p, RNG& rng)
{
  size_t offset=metapopulation_idx*nodes_per_metapopulation;

  std::vector<size_t> connectable;
  for (size_t node_idx=offset;
      node_idx<offset+nodes_per_metapopulation;
      ++node_idx)
  {
    auto degree=out_degree(node_idx, g);
    boost::random::binomial_distribution<> binomial(degree, binomial_p);
    auto connection_cnt=binomial(rng);
    for (size_t connect_idx=0; connect_idx<connection_cnt; ++connect_idx)
    {
      connectable.push_back(node_idx);
    }
  }

  return connectable;
}




/*! Construct a sequence of erdos-renyi metapopulations connected
 *  with binomial samples on the edge degree of each node.
 */
template<typename RNG>
ContactGraph serial_metapopulation(size_t nodes_per_metapopulation,
    size_t metapopulation_cnt, double edge_fraction, double binomial_p,
    RNG& rng)
{
  ContactGraph g(nodes_per_metapopulation*metapopulation_cnt);

  // Make multiple independent metapopulations, each from its own
  // erdos_renyi distribution.
  using ERGen=boost::erdos_renyi_iterator<RNG, ContactGraph>;
  for (size_t make_meta_idx=0; make_meta_idx<metapopulation_cnt; ++make_meta_idx)
  {
    size_t offset=make_meta_idx*nodes_per_metapopulation;

    ERGen erdos_renyi_iter(rng, nodes_per_metapopulation,
        edge_fraction);
    size_t edge_cnt=0;
    for (; erdos_renyi_iter!=ERGen(); ++erdos_renyi_iter)
    {
      size_t source_vertex, target_vertex;
      std::tie(source_vertex, target_vertex)=*erdos_renyi_iter;
      // remove self-loops here.
      if (source_vertex!=target_vertex)
      {
        add_edge(source_vertex+offset, target_vertex+offset, g);
        ++edge_cnt;
      }
    }
    BOOST_LOG_TRIVIAL(info)<<"metapop "<<make_meta_idx<<" has "
      <<edge_cnt<<" edges";
  }

  // Connect the distinct metapopulations.
  for (size_t left_idx=0; left_idx<metapopulation_cnt-1; ++left_idx)
  {
    std::vector<size_t> left_connectable=
      available_to_connect(g, left_idx, nodes_per_metapopulation,
      binomial_p, rng);
    std::vector<size_t> right_connectable=
      available_to_connect(g, left_idx+1, nodes_per_metapopulation,
      binomial_p, rng);

    size_t connection_cnt;
    if (left_connectable.size() > right_connectable.size())
    {
        std::shuffle(left_connectable.begin(), left_connectable.end(), rng);
        connection_cnt=right_connectable.size();
    }
    else
    {
        std::shuffle(right_connectable.begin(), right_connectable.end(), rng);
        connection_cnt=left_connectable.size();
    }

    for (size_t cidx=0; cidx<connection_cnt; ++cidx)
    {
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
  inline friend std::ostream&
  operator<<(std::ostream& os, const IndividualToken& it)
  {
    return os << 't';
  }
};



struct SIRPlace
{
  size_t disease;
  size_t individual;
  size_t metapop;

  SIRPlace()=default;
  SIRPlace(size_t d, size_t i, size_t m)
  : disease(d), individual(i), metapop(m)
  {}

  friend inline
  bool operator<(const SIRPlace& a, const SIRPlace& b)
  {
    return lazy_less(a.disease, b.disease, a.individual,
      b.individual, a.metapop, b.metapop);
  }


  friend inline
  bool operator==(const SIRPlace& a, const SIRPlace& b)
  {
    return (a.disease==b.disease)&& (a.individual==b.individual)
        && (a.metapop==b.metapop);
  }


  friend inline
  std::ostream&
  operator<<(std::ostream& os, const SIRPlace& cp)
  {
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
  bool operator<(const SIRTKey& a, const SIRTKey& b)
  {
    return lazy_less(a.ind1, b.ind1, a.ind2, b.ind2);
  }

  friend inline
  bool operator==(const SIRTKey& a, const SIRTKey& b)
  {
    return (a.ind1==b.ind1) && (a.ind2==b.ind2);
  }

  friend inline
  std::ostream&
  operator<<(std::ostream& os, const SIRTKey& cp)
  {
    return os << '(' << cp.ind1 << ", " << cp.ind2 << ')';
  }
};



using Local=LocalMarking<Uncolored<IndividualToken>>;
// Marking of the net.
using Mark=Marking<size_t, Uncolored<IndividualToken>>;
// State of the continuous dynamical system.
class WithParams
{
  std::map<int,double> params;
};

using SIRState=GSPNState<Mark,WithParams>;

using SIRGSPN=ExplicitTransitions<SIRPlace, SIRTKey, Local, RNG, WithParams>;


using Dist=TransitionDistribution<RNG>;
using ExpDist=ExponentialDistribution<RNG>;

using SIRTransition=ExplicitTransition<Local,RNG,WithParams>;


// Now make specific transitions.
class InfectNeighbor : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  enabled(const WithParams& s, const Local& lm,
    double te, double t0) const override
  {
    if (lm.template input_tokens_sufficient<0>())
    {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    }
    else
    {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void fire(WithParams& s, Local& lm,
      RNG& rng) const override
  {
    BOOST_LOG_TRIVIAL(debug) << "Fire infection " << lm;
    lm.template transfer_by_stochiometric_coefficient<0>(rng);
  }

};





// Now make specific transitions.
class Recover : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  enabled(const WithParams& s, const Local& lm,
      double te, double t0) const override
  {
    if (lm.template input_tokens_sufficient<0>())
    {
      return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0, te))};
    }
    else
    {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void fire(WithParams& s, Local& lm,
      RNG& rng) const override
  {
    BOOST_LOG_TRIVIAL(debug) << "Fire recovery "<< lm;
    lm.template transfer_by_stochiometric_coefficient<0>(rng);
  }
};




SIRGSPN
build_system(size_t nodes_per_metapopulation,
    size_t metapopulation_cnt, double edge_fraction, double binomial_p,
    RNG& rng)
{
  auto contact_graph=serial_metapopulation(nodes_per_metapopulation,
    metapopulation_cnt, edge_fraction, binomial_p, rng);

  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  enum { s, i, r };

  for (size_t meta_idx=0; meta_idx<metapopulation_cnt; ++meta_idx)
  {
    for (size_t ind_idx=0; ind_idx<nodes_per_metapopulation; ++ind_idx)
    {
      auto sp=SIRPlace{s, ind_idx, meta_idx};
      bg.add_place(sp);
      auto ip=SIRPlace{i, ind_idx, meta_idx};
      bg.add_place(ip);
      auto rp=SIRPlace{r, ind_idx, meta_idx};
      bg.add_place(rp);

      bg.add_transition({ip, rp},
        {Edge{ip, -1}, Edge{rp, -1}},
        std::unique_ptr<SIRTransition>(new Recover()));
    }
  }

  auto index_to_place=[&](size_t index)->SIRPlace
  {
    size_t mpop=index / nodes_per_metapopulation;
    size_t individual=index % nodes_per_metapopulation;
    return SIRPlace(s, individual, mpop);
  };

  std::set<std::array<size_t,2>> seen;
  for (auto ce=edges(contact_graph); ce.first!=ce.second; ++ce.first)
  {
    auto left_index=source(*ce.first, contact_graph);
    auto lefts=index_to_place(left_index);
    auto right_index=target(*ce.first, contact_graph);
    auto rights=index_to_place(right_index);

    std::array<size_t,2> new_transition{
      std::min(left_index, right_index), std::max(left_index, right_index)};
    if (seen.find(new_transition)==seen.end())
    {
      SIRPlace lefti{i, lefts.individual, lefts.metapop};
      SIRPlace righti{i, rights.individual, rights.metapop};

      bg.add_transition({lefti, rights},
        {Edge{lefti, -1}, Edge{rights, -1}, Edge{lefti, 1}, Edge{righti, 1}},
        std::unique_ptr<SIRTransition>(new InfectNeighbor()));

      bg.add_transition({righti, lefts},
        {Edge{righti, -1}, Edge{lefts, -1}, Edge{righti, 1}, Edge{lefti, 1}},
        std::unique_ptr<SIRTransition>(new InfectNeighbor()));

      seen.insert(new_transition);
    }
  }

  return std::move(bg.build());
}





template<typename GSPN>
class SIROutputFunction
{
public:
  typedef void result_type;

private:
  size_t _pop_cnt;
  size_t _ind_cnt;
  size_t _threshold;
  std::vector<size_t> _infected_per_pop;
  std::vector<double> _first_passage_time;
  std::vector<bool> _passed;

public:
  SIROutputFunction(size_t population_cnt,
    size_t individuals_per_metapopulation, size_t threshold_for_passage)
  : _pop_cnt(population_cnt), _ind_cnt(individuals_per_metapopulation),
    _threshold(threshold_for_passage),
    _first_passage_time(population_cnt, 0), _passed(population_cnt, false),
    _infected_per_pop(population_cnt, 0)
  {

  }



  result_type operator()(const GSPN& gspn, const SIRState& state)
  {
    auto& modified=state.marking.modified();
    for (auto place_idx : modified)
    {
      SIRPlace p=gspn.vertex_place(place_idx);
      if (p.disease==1)
      {
        bool filled=(length<0>(state.marking, place_idx)>0);
        if (filled)
        {
          _infected_per_pop[p.metapop]+=1;
          if (_infected_per_pop[p.metapop]==_threshold)
          {
            _first_passage_time[p.metapop]=state.current_time();
            _passed[p.metapop]=true;
            BOOST_LOG_TRIVIAL(info)<<"Population "<<p.metapop
              <<" at time "<<state.current_time();
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
  size_t individuals_per_metapopulation=10000;
  size_t metapopulation_cnt=20;
  size_t rand_seed=1;
  std::string log_level;

  desc.add_options()
    ("help", "show help message")
    ("size,s",
      po::value<size_t>(&individuals_per_metapopulation)->default_value(1000),
      "size of each partially-mixed population")
    ("seed,r",
      po::value<size_t>(&rand_seed)->default_value(1),
      "seed for random number generator")
    ("count,c",
      po::value<size_t>(&metapopulation_cnt)->default_value(20),
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

  afidd::log_init(log_level);
  RNG rng(rand_seed);


  auto gspn=build_system(individuals_per_metapopulation, metapopulation_cnt,
    poisson_constant/individuals_per_metapopulation, binomial_to_neighbors, rng);

  size_t individual_cnt=metapopulation_cnt*individuals_per_metapopulation;
  BOOST_LOG_TRIVIAL(debug)<<"Initializing "<<individual_cnt<< " individuals "
      "as susceptible.";

  SIRState state;
  for (size_t meta_idx=0; meta_idx<metapopulation_cnt; ++meta_idx)
  {
    for (size_t ind_idx=0; ind_idx<individuals_per_metapopulation; ++ind_idx)
    {
      auto vert=gspn.place_vertex(SIRPlace{0, ind_idx, meta_idx});
      add<0>(state.marking, vert, IndividualToken{});
    }
  }

  using Markov=PartialCoreMatrix<SIRGSPN, SIRState, RNG>;
  Markov system(gspn, state);

  SIROutputFunction<SIRGSPN> output(
    metapopulation_cnt, individuals_per_metapopulation,
    std::lround(0.1*individuals_per_metapopulation));

  // The initial input string moves a token from susceptible to infected.
  auto first_case=afidd::smv::uniform_index(rng, individuals_per_metapopulation);
  auto susceptible=gspn.place_vertex(SIRPlace{0, first_case, 0});
  auto infected=gspn.place_vertex(SIRPlace{1, first_case, 0});

  auto input_string=[&susceptible, &infected](SIRState& state)->void {
    move<0,0>(state.marking, susceptible, infected, 1);
  };
  auto next=propagate_competing_processes(system, input_string, rng);

  output(gspn, state);

  size_t transition_cnt=0;
  auto nothing=[](SIRState&)->void {};
  for ( ;
    std::get<1>(next)<std::numeric_limits<double>::max();
    next=propagate_competing_processes(system, nothing, rng))
  {
    BOOST_LOG_TRIVIAL(debug) << "trans " << std::get<0>(next) << " time " <<
        std::get<1>(next);
    ++transition_cnt;
    output(gspn, state);
  }
  BOOST_LOG_TRIVIAL(info)<< transition_cnt << " transitions";
  return 0;
}

