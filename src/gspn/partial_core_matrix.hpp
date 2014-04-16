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
#include "continuous_dynamics.hpp"

namespace afidd
{
namespace smv
{

template<typename GSPN, typename State, typename RNG>
class PartialCoreMatrix
{
 public:
  // Re-advertise the transition key.
  typedef typename GSPN::TransitionKey TransitionKey;
  typedef ContinuousPropagator<TransitionKey,RNG> Propagator;
  using PropagatorVector=std::vector<Propagator*>;
  typedef GSPN PetriNet;

  PartialCoreMatrix(GSPN& gspn, PropagatorVector pv)
  : gspn_(gspn), propagator_{pv} {}

  void set_state(State* s) {
    state_=s;
  }

  void MakeCurrent(RNG& rng) {
    if (state_->marking.Modified().size()==0) return;
    // Check all neighbors of a place to see if they were enabled.
    auto lm=state_->marking.GetLocalMarking();

    NeighborsOfPlaces(gspn_, state_->marking.Modified(),
      [&] (TransitionKey neighbor_id) {
        // Was this transition enabled? When?
        double enabling_time=0.0;
        Propagator* previous_propagator=nullptr;
        for (const auto& enable_prop : propagator_) {
          double previously_enabled=false;
          std::tie(previously_enabled, enabling_time)
              =enable_prop->Enabled(neighbor_id);
          if (previously_enabled) {
            previous_propagator=enable_prop;
            break;
          }
        }
        if (previous_propagator==nullptr) enabling_time=state_->CurrentTime();

        // Set up the local marking.
        auto neighboring_places=
            NeighborsOfTransition(gspn_, neighbor_id);
        state_->marking.InitLocal(lm, neighboring_places);

        bool isEnabled=false;
        std::unique_ptr<TransitionDistribution<RNG>> dist;
        std::tie(isEnabled, dist)=
            Enabled(gspn_, neighbor_id, state_->user, lm,
            enabling_time, state_->CurrentTime());

        if (isEnabled) {
          Propagator* appropriate=nullptr;
          for (const auto& prop_ptr : propagator_) {
            if (prop_ptr->Include(*dist)) {
              appropriate=prop_ptr;
            }
          }
          BOOST_ASSERT_MSG(appropriate!=nullptr, "No propagator willing to "
              "accept this distribution");
          // Even if it was already enabled, take the new distribution
          // in case it has changed.
          if (dist!=nullptr) {
            bool was_enabled=previous_propagator!=nullptr;
            if (was_enabled) {
              if (previous_propagator==appropriate) {
                appropriate->Enable(neighbor_id, dist, state_->CurrentTime(),
                    was_enabled, rng);
              } else {
                previous_propagator->Disable(neighbor_id, state_->CurrentTime());
                appropriate->Enable(neighbor_id, dist, state_->CurrentTime(),
                    was_enabled, rng);
              }
            } else {
              appropriate->Enable(neighbor_id, dist, state_->CurrentTime(),
                  was_enabled, rng);
            }
          } else {
            BOOST_ASSERT_MSG(previous_propagator!=nullptr, "Transition didn't "
                "return a distribution, so it thinks it was enabled, but it "
                "isn't listed as enabled in any propagator");
          }
        } else if (!isEnabled && previous_propagator!=nullptr) {
          previous_propagator->Disable(neighbor_id, state_->CurrentTime());

        } else {
          ; // not enabled, not becoming enabled.
        }
      });
    BOOST_LOG_TRIVIAL(trace) << "Marking modified cnt: "<<
        state_->marking.Modified().size();
    state_->marking.Clear();
  }


  void Trigger(TransitionKey trans_id, double when, RNG& rng) {
    auto neighboring_places=NeighborsOfTransition(gspn_, trans_id);

    auto lm=state_->marking.GetLocalMarking();
    state_->marking.InitLocal(lm, neighboring_places);
    Fire(gspn_, trans_id, state_->user, lm, rng);
    state_->marking.ReadLocal(lm, neighboring_places);

    BOOST_LOG_TRIVIAL(trace) << "Fire "<<trans_id <<" neighbors: "<<
        neighboring_places.size() << " modifies "
        << state_->marking.Modified().size() << " places.";

    auto current_time=state_->SetTime(when);

    bool enabled=false;
    double previous_when;
    for (auto& prop_ptr : propagator_) {
      std::tie(enabled, previous_when)=prop_ptr->Enabled(trans_id);
      if (enabled) {
        prop_ptr->Fire(trans_id, state_->CurrentTime(), rng);
        break;
      }
    }
    BOOST_ASSERT_MSG(enabled, "The transition that fired wasn't enabled?");
  }

  PropagatorVector& Propagators() { return propagator_; }

 private:
  GSPN& gspn_;
  State* state_;
  PropagatorVector propagator_;
};

} // smv
} // afidd

#endif // _PARTIAL_CORE_MATRIX_H_
