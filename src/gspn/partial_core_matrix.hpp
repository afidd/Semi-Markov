#ifndef _PARTIAL_CORE_MATRIX_H_
#define _PARTIAL_CORE_MATRIX_H_ 1

#include <memory>
#include <map>
#include "logging.hpp"
#include "distributions.hpp"

namespace afidd
{
namespace smv
{

template<typename GSPN, typename State, typename RNG>
class PartialCoreMatrix
{
  GSPN& _gspn;
  using Marking=typename State::Marking;
  State& _state;
  using Dist=TransitionDistribution<RNG>;
  std::map<typename GSPN::TransitionKey,std::unique_ptr<Dist>> _distributions;

public:
  typedef GSPN PetriNet;
  // Re-advertise the transition key.
  typedef typename GSPN::TransitionKey TransitionKey;

  PartialCoreMatrix(GSPN& gspn, State& s)
  : _gspn(gspn), _state(s)
  {}



  template<typename FUNCTOR>
  void state_machine_token(const FUNCTOR& token)
  {
    token(_state);
  }



  template<typename FUNCTOR>
  void transitions(const FUNCTOR& eval)
  { 
    if (_state.marking.Modified().size()>0)
    {
      // Check all neighbors of a place to see if they were enabled.
      auto lm=_state.marking.GetLocalMarking();

      neighbors_of_places(_gspn, _state.marking.Modified(),
        [&] (TransitionKey neighbor_id)
        {
          // Was this transition enabled? When?
          auto previous_distribution=_distributions.find(neighbor_id);
          double enabling_time;
          bool previously_enabled=previous_distribution!=_distributions.end();
          if (previously_enabled)
          {
            enabling_time=previous_distribution->second->EnablingTime();
          }
          else
          {
            enabling_time=_state.CurrentTime();
          }

          // Set up the local marking.
          auto neighboring_places=
              neighbors_of_transition(_gspn, neighbor_id);
          _state.marking.InitLocal(lm, neighboring_places);

          bool isEnabled=false;
          std::unique_ptr<TransitionDistribution<RNG>> dist;
          std::tie(isEnabled, dist)=enabled(_gspn, neighbor_id, _state.user, lm,
              enabling_time, _state.CurrentTime());

          if (isEnabled)
          {
            // Even if it was already enabled, take the new distribution
            // in case it has changed.
            if (dist!=nullptr)
            {
              _distributions.emplace(neighbor_id, std::move(dist));
            }
            // else it's OK if they return nullptr. Use old distribution.
          }
          else if (!isEnabled && previously_enabled)
          {
            _distributions.erase(neighbor_id);
          }
          else
          {
            ; // not enabled, not becoming enabled.
          }
        });
      BOOST_LOG_TRIVIAL(trace) << "Marking modified cnt: "<<
          _state.marking.Modified().size() << " enabled " <<
          _distributions.size();
      _state.marking.Clear();
    }

    auto begin=_distributions.begin();
    for (; begin!=_distributions.end(); ++begin)
    {
      TransitionKey trans_id=begin->first;
      eval(begin->second, trans_id, _state.CurrentTime());
    }
  }


  void trigger(TransitionKey trans_id, double when, RNG& rng)
  {
    auto neighboring_places=neighbors_of_transition(_gspn, trans_id);

    BOOST_LOG_TRIVIAL(trace)<<"Fire "<<trans_id<<" neighbors: "<<
      neighboring_places.size();
    auto lm=_state.marking.GetLocalMarking();
    _state.marking.InitLocal(lm, neighboring_places);
    fire(_gspn, trans_id, _state.user, lm, rng);
    _state.marking.ReadLocal(lm, neighboring_places);

    BOOST_LOG_TRIVIAL(trace) << "Fire "<<trans_id << " modifies "
      << _state.marking.Modified().size() << " places.";

    auto current_time=_state.AddTime(when);
    _distributions.erase(trans_id);
  }
};

} // smv
} // afidd

#endif // _PARTIAL_CORE_MATRIX_H_
