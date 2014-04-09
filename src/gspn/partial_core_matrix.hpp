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
  GSPN& gspn_;
  using Marking=typename State::Marking;
  State& state_;
  using Dist=TransitionDistribution<RNG>;
  std::map<typename GSPN::TransitionKey,std::unique_ptr<Dist>> distributions_;

public:
  typedef GSPN PetriNet;
  // Re-advertise the transition key.
  typedef typename GSPN::TransitionKey TransitionKey;

  PartialCoreMatrix(GSPN& gspn, State& s)
  : gspn_(gspn), state_(s)
  {}


  template<typename FUNCTOR>
  void StateMachineToken(const FUNCTOR& token)
  {
    token(state_);
  }


  template<typename FUNCTOR>
  void Transitions(const FUNCTOR& eval)
  { 
    if (state_.marking.Modified().size()>0)
    {
      // Check all neighbors of a place to see if they were enabled.
      auto lm=state_.marking.GetLocalMarking();

      NeighborsOfPlaces(gspn_, state_.marking.Modified(),
        [&] (TransitionKey neighbor_id)
        {
          // Was this transition enabled? When?
          auto previous_distribution=distributions_.find(neighbor_id);
          double enabling_time;
          bool previously_enabled=previous_distribution!=distributions_.end();
          if (previously_enabled)
          {
            enabling_time=previous_distribution->second->EnablingTime();
          }
          else
          {
            enabling_time=state_.CurrentTime();
          }

          // Set up the local marking.
          auto neighboring_places=
              NeighborsOfTransition(gspn_, neighbor_id);
          state_.marking.InitLocal(lm, neighboring_places);

          bool isEnabled=false;
          std::unique_ptr<TransitionDistribution<RNG>> dist;
          std::tie(isEnabled, dist)=
              Enabled(gspn_, neighbor_id, state_.user, lm,
              enabling_time, state_.CurrentTime());

          if (isEnabled)
          {
            // Even if it was already enabled, take the new distribution
            // in case it has changed.
            if (dist!=nullptr)
            {
              distributions_.emplace(neighbor_id, std::move(dist));
            }
            // else it's OK if they return nullptr. Use old distribution.
          }
          else if (!isEnabled && previously_enabled)
          {
            distributions_.erase(neighbor_id);
          }
          else
          {
            ; // not enabled, not becoming enabled.
          }
        });
      BOOST_LOG_TRIVIAL(trace) << "Marking modified cnt: "<<
          state_.marking.Modified().size() << " enabled " <<
          distributions_.size();
      state_.marking.Clear();
    }

    auto begin=distributions_.begin();
    for (; begin!=distributions_.end(); ++begin)
    {
      TransitionKey trans_id=begin->first;
      eval(begin->second, trans_id, state_.CurrentTime());
    }
  }


  void Trigger(TransitionKey trans_id, double when, RNG& rng)
  {
    auto neighboring_places=NeighborsOfTransition(gspn_, trans_id);

    BOOST_LOG_TRIVIAL(trace)<<"Fire "<<trans_id<<" neighbors: "<<
      neighboring_places.size();
    auto lm=state_.marking.GetLocalMarking();
    state_.marking.InitLocal(lm, neighboring_places);
    Fire(gspn_, trans_id, state_.user, lm, rng);
    state_.marking.ReadLocal(lm, neighboring_places);

    BOOST_LOG_TRIVIAL(trace) << "Fire "<<trans_id << " modifies "
      << state_.marking.Modified().size() << " places.";

    auto current_time=state_.AddTime(when);
    distributions_.erase(trans_id);
  }
};

} // smv
} // afidd

#endif // _PARTIAL_CORE_MATRIX_H_
