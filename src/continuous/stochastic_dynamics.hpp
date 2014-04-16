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
#ifndef _STOCHASTIC_DYNAMICS_H_
#define _STOCHASTIC_DYNAMICS_H_ 1

#include <tuple>
#include <limits>
#include "partial_core_matrix.hpp"

namespace afidd
{
namespace smv
{
template<typename GSPN,typename DynState,typename RandGen>
class StochasticDynamics
{
  using Kernel=PartialCoreMatrix<GSPN,DynState,RandGen>;
  using PropagatorVector=typename Kernel::PropagatorVector;
 public:
  // Re-advertise the state.
  using State=DynState;
  StochasticDynamics(GSPN& gspn, PropagatorVector pv)
  : core_matrix_(gspn, pv) {}

  void Initialize(State* state, RandGen *rng) {
    core_matrix_.set_state(state);
    rng_=rng;
  }

  bool operator()(State& state) {
    core_matrix_.MakeCurrent(*rng_);
    using TransitionKey=typename Kernel::TransitionKey;
    std::tuple<TransitionKey,double> least{TransitionKey{},
        std::numeric_limits<double>::infinity()};
    for (auto& propagator : core_matrix_.Propagators()) {
      auto next=propagator->Next(state.CurrentTime(), *rng_);
      if (std::get<1>(next)<std::get<1>(least)) {
        least=next;
      }
    }
    bool running=true;
    if (std::get<1>(least)<(std::numeric_limits<double>::max)()) {
      core_matrix_.Trigger(std::get<0>(least), std::get<1>(least), *rng_);
    } else {
      running=false;
    }
    state.last_transition=std::get<0>(least);
    return running;
  }
 private:
  Kernel core_matrix_;
  RandGen *rng_;
};

}
}

#endif // _STOCHASTIC_DYNAMICS_H_ //
