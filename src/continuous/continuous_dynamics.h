#ifndef _CONTINUOUS_DYNAMICS_H_
#define _CONTINUOUS_DYNAMICS_H_ 1

namespace afidd
{
namespace smv
{
template<typename EmbeddedMarkov, typename T, typename RNG>
std::tuple<trans_t<typename EmbeddedMarkov::PetriNet>, double>
delta(EmbeddedMarkov& system, T& token, RNG& rng)
{
  system.state_machine_token(token);
  using Transition=trans_t<typename EmbeddedMarkov::PetriNet>;

  auto least=std::make_tuple(Transition{}, std::numeric_limits<double>::infinity());

  using DistPtr=std::unique_ptr<TransitionDistribution<RNG>>;

  system.transitions(
    [&least, &rng] (std::unique_ptr<TransitionDistribution<RNG>>& distribution,
          Transition trans_id, double enabling_time, double now) {
      auto trial_time=distribution->sample(enabling_time, now, rng);
      if (trial_time < std::get<1>(least))
      {
        std::get<0>(least)=trans_id;
        std::get<1>(least)=trial_time;
      }
    });

  if (std::get<1>(least)<std::numeric_limits<double>::infinity())
  {
    system.fire(std::get<0>(least), std::get<1>(least), rng);
  }
  return least;
}

} // smv
} // afidd


#endif // _CONTINUOUS_DYNAMICS_H_
