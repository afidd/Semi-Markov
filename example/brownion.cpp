/*! An example that doesn't use explicit transitions.
 *  The two-state Brownion is a particle on a 2d rectangular lattice
 *  that has two sets of transition distributions depending on its state.
 *  The lattice is infinite in both directions.
 */

#include "stochnet.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/program_options.hpp"
#include "partial_core_matrix.hpp"
#include "continuous_dynamics.hpp"
#include "logging.hpp"
#include "brownion_model.hpp"



namespace smv=afidd::smv;


int main(int argc, char *argv[])
{
  size_t iteration_cnt=100;

  namespace po=boost::program_options;
  po::options_description desc("Two-state Brownion.");
  size_t rand_seed=1;
  double beta=1.0;
  double gamma=1.0;
  std::string log_level;

  desc.add_options()
    ("help", "show help message")
    ("size,s",
      po::value<size_t>(&iteration_cnt)->default_value(100),
      "number of steps to take")
    ("seed,r",
      po::value<size_t>(&rand_seed)->default_value(1),
      "seed for random number generator")
    ("beta",
      po::value<double>(&beta)->default_value(1.0),
      "parameter")
    ("gamma",
      po::value<double>(&gamma)->default_value(1.0),
      "parameter")
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

  RandGen rng(rand_seed);
  BrownionGSPN gspn;
  BrownionState state;

  using Markov=PartialCoreMatrix<BrownionGSPN,BrownionState,RandGen>;
  Markov system(gspn, state);

  auto input_string=[](BrownionState& state)->void {
    Add<0>(state.marking, PlaceType{0,0,0}, IndividualToken());
  };
  auto nothing=[](BrownionState&)->void {};

  auto next=propagate_competing_processes(system, input_string, rng);
  for (size_t iteration_idx=0; iteration_idx<iteration_cnt; ++iteration_idx)
  {
    BOOST_LOG_TRIVIAL(info) << "trans " << std::get<0>(next) << " time " <<
        std::get<1>(next);
    next=propagate_competing_processes(system, nothing, rng);
  }

  return 0;
}
