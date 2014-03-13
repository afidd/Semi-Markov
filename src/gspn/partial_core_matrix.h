#ifndef _PARTIAL_CORE_MATRIX_H_
#define _PARTIAL_CORE_MATRIX_H_ 1

#include <memory>
#include <map>
#include "logging.h"

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
  std::map<trans_t<GSPN>,std::unique_ptr<Dist>> _distributions;

public:
  typedef GSPN PetriNet;
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
    if (_state.marking.modified().size()>0)
    {
      // Check all neighbors of a place to see if they were enabled.
      neighbors_of_places(_gspn, _state.marking.modified(),
        [&] (trans_t<GSPN> neighbor_id)
        {
          auto neighboring_places=
              neighbors_of_transition(_gspn, neighbor_id);
          LocalMarking<Marking> lm(_state.marking, neighboring_places);

          bool isEnabled;
          std::unique_ptr<TransitionDistribution<RNG>> dist;
          std::tie(isEnabled, dist)=enabled(_gspn, neighbor_id, _state, lm);
          auto previously_enabled=
            (_state.enabling_time.find(neighbor_id)!=_state.enabling_time.end());
          if (isEnabled && !previously_enabled)
          {
            _state.enabling_time.emplace(neighbor_id, _state.current_time());
            _distributions.emplace(neighbor_id, std::move(dist));
          }
          else if (!isEnabled && previously_enabled)
          {
            _state.enabling_time.erase(neighbor_id);
            _distributions.erase(neighbor_id);
          }
          else
          {
            ; // No change to this distribution.
          }
        });
      BOOST_LOG_TRIVIAL(trace) << "Marking modified cnt: "<<
          _state.marking.modified().size() << " enabled " <<
          _distributions.size();
      _state.marking.clear();
    }

    auto begin=_distributions.begin();
    for (; begin!=_distributions.end(); ++begin)
    {
      trans_t<GSPN> trans_id=begin->first;
      eval(begin->second, trans_id, _state.enabling_time[trans_id],
        _state.current_time());
    }
  }


  void fire(trans_t<GSPN> trans_id, double when, RNG& rng)
  {
    auto neighboring_places=neighbors_of_transition(_gspn, trans_id);
    LocalMarking<Marking> lm(_state.marking, neighboring_places);
    afidd::smv::fire(_gspn, trans_id, _state, lm, rng);
    BOOST_LOG_TRIVIAL(trace) << "fire "<<trans_id << " modifies "
      << _state.marking.modified().size() << " places.";

    auto current_time=_state.add_time(when);
    _state.enabling_time.erase(trans_id);
    _distributions.erase(trans_id);
  }
};

} // smv
} // afidd

#endif // _PARTIAL_CORE_MATRIX_H_
