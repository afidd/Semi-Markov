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
#ifndef _CONTINUOUS_DYNAMICS_H_
#define _CONTINUOUS_DYNAMICS_H_ 1

#include <tuple>
#include <limits>
#include <memory>
#include "distributions.hpp"


namespace afidd
{
namespace smv
{
template<typename PartialCore, typename T, typename RNG>
std::tuple<typename PartialCore::TransitionKey, double>
PropagateCompetingProcesses(PartialCore& system, T& token, RNG& rng)
{
  system.StateMachineToken(token);
  using Transition=typename PartialCore::TransitionKey;

  auto least=std::make_tuple(Transition{},
      std::numeric_limits<double>::infinity());

  using DistPtr=std::unique_ptr<TransitionDistribution<RNG>>;

  system.Transitions(
    [&least, &rng] (std::unique_ptr<TransitionDistribution<RNG>>& distribution,
          Transition trans_id, double now) {
      auto trial_time=distribution->Sample(now, rng);
      if (trial_time < std::get<1>(least))
      {
        std::get<0>(least)=trans_id;
        std::get<1>(least)=trial_time;
      }
    });

  if (std::get<1>(least)<std::numeric_limits<double>::infinity()) {
    system.Trigger(std::get<0>(least), std::get<1>(least), rng);
  }
  return least;
}

} // smv
} // afidd


#endif // _CONTINUOUS_DYNAMICS_H_
