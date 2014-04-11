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
#include "boost/utility.hpp"
#include "distributions.hpp"


namespace afidd
{
namespace smv
{


template<typename RNG>
struct IncludeDistributions
{
  virtual bool Include(const TransitionDistribution<RNG>& distribution) const {
    return true;
  }
};


/*! Store a map of every enabled transition.
 *  The interface is meant for the PartialCoreMatrix.
 */
template<typename TransitionKey,typename RNG>
class TransitionsStoreEnabled : boost::noncopyable
{
 public:
  TransitionsStoreEnabled() {}
  ~TransitionsStoreEnabled() {}

  std::tuple<bool,double> Enabled(const TransitionKey& tkey) const {
    auto previous_distribution=distributions_.find(tkey);
    bool previously_enabled=previous_distribution!=distributions_.end();
    if (previously_enabled) {
      double enabling_time=previous_distribution->second->EnablingTime();
      return std::make_tuple(true, enabling_time);
    } else {
      return std::make_tuple(false, 0.0);
    }
  }

  void Enable(const TransitionKey& tkey,
      std::unique_ptr<TransitionDistribution<RNG>>& distribution,
      bool previously_enabled) {
    distributions_.emplace(tkey, std::move(distribution));
  }

  void Disable(const TransitionKey& tkey) {
    distributions_.erase(tkey);
  }

  template<typename EachTransition>
  void All(const EachTransition& eval) const {
    auto begin=distributions_.begin();
    for (; begin!=distributions_.end(); ++begin) {
      TransitionKey trans_id=begin->first;
      eval(begin->second, trans_id);
    }
  }

 private:
  using Dist=TransitionDistribution<RNG>;
  std::map<TransitionKey,std::unique_ptr<Dist>> distributions_;
};


template<typename TransitionKey,typename RNG>
class ContinuousPropagator
{
 public:
  virtual bool Include(const TransitionDistribution<RNG>& distribution) const=0;
  virtual std::tuple<TransitionKey, double> Next(double now, RNG& rng) const=0;
  virtual std::tuple<bool,double> Enabled(const TransitionKey& tkey) const=0;
  virtual void Enable(const TransitionKey& tkey,
    std::unique_ptr<TransitionDistribution<RNG>>& distribution,
    bool previously_enabled)=0;
  virtual void Disable(const TransitionKey& tkey)=0;
};


template<typename TransitionKey,typename RNG>
class PropagateCompetingProcesses
: public ContinuousPropagator<TransitionKey,RNG>
{
 public:
  PropagateCompetingProcesses() : include_{new IncludeDistributions<RNG>{}} {}

  bool Include(const TransitionDistribution<RNG>& distribution) const {
    return include_->Include(distribution);
  }

  std::tuple<TransitionKey, double> Next(double now, RNG& rng) const {
    auto least=std::make_tuple(TransitionKey{},
        std::numeric_limits<double>::infinity());

    using DistPtr=std::unique_ptr<TransitionDistribution<RNG>>;

    transitions_.All(
      [&least, &rng, &now] (const DistPtr& distribution,
          TransitionKey trans_id)->void {
        auto trial_time=distribution->Sample(now, rng);
        if (trial_time < std::get<1>(least)) {
          std::get<0>(least)=trans_id;
          std::get<1>(least)=trial_time;
        }
      });

    return least;
  }

  std::tuple<bool,double> Enabled(const TransitionKey &tkey) const {
    return transitions_.Enabled(tkey);
  }

  void Enable(const TransitionKey &tkey,
      std::unique_ptr<TransitionDistribution<RNG>> &distribution,
      bool previously_enabled) {
    transitions_.Enable(tkey, distribution, previously_enabled);
  }

  void Disable(const TransitionKey& tkey) {
    transitions_.Disable(tkey);
  }

 private:
  TransitionsStoreEnabled<TransitionKey,RNG> transitions_;
  std::unique_ptr<IncludeDistributions<RNG>> include_;
};

} // smv
} // afidd


#endif // _CONTINUOUS_DYNAMICS_H_
