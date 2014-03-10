/*! An example that doesn't use explicit transitions.
 *  The two-state Brownion is a particle on a 2d rectangular lattice
 *  that has two sets of transition distributions depending on its state.
 *  The lattice is infinite in both directions.
 */

#include "stochnet.h"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/mpl/vector.hpp"
#include "brownion_model.h"
#include "embedded_markov.h"
#include "continuous_dynamics.h"
#include "logging.h"



namespace smv=afidd::smv;


int main(int argc, char *argv[])
{
  size_t iteration_cnt=100;
  size_t random_seed=1;
  afidd::log_init("debug");

  RandGen rng(random_seed);
  BrownionGraph graph;
  Transitions transitions;
  BrownionState state;

  using Markov=EmbeddedMarkov<BrownionGraph,Transitions,BrownionState,RandGen>;
  Markov system(graph, transitions, state);

  auto input_string=[](BrownionState& state)->void {
    add<0>(state.marking, PlaceType{0,0,0}, IndividualToken());
  };
  auto nothing=[](BrownionState&)->void {};

  auto next=delta(system, input_string, rng);
  for (size_t iteration_idx=0; iteration_idx<iteration_cnt; ++iteration_idx)
  {
    std::get<1>(next)<std::numeric_limits<double>::max();
    BOOST_LOG_TRIVIAL(debug) << "trans " << std::get<0>(next) << " time " <<
        std::get<1>(next);
    next=delta(system, nothing, rng);
  }

  return 0;
}
