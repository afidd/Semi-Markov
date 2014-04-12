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
#include "boost/heap/fibonacci_heap.hpp"
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


template<typename RNG>
struct BoundedHazardDistributions : public IncludeDistributions<RNG>
{
  virtual bool Include(const TransitionDistribution<RNG>& distribution) const {
    return distribution.BoundedHazard();
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
    double when,
    bool previously_enabled,
    RNG& rng)=0;
  virtual void Disable(const TransitionKey& tkey, double when)=0;
  virtual void Fire(const TransitionKey& tkey, double when, RNG& rng)=0;
};


template<typename TransitionKey,typename RNG>
class PropagateCompetingProcesses
: public ContinuousPropagator<TransitionKey,RNG>
{
 public:
  PropagateCompetingProcesses() : include_{new IncludeDistributions<RNG>{}} {}

  virtual bool Include(const TransitionDistribution<RNG>& distribution
      ) const override {
    return include_->Include(distribution);
  }

  virtual std::tuple<TransitionKey, double> Next(double now, RNG& rng
      ) const override {
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

  virtual std::tuple<bool,double> Enabled(const TransitionKey &tkey
      ) const override {
    return transitions_.Enabled(tkey);
  }

  virtual void Enable(const TransitionKey &tkey,
      std::unique_ptr<TransitionDistribution<RNG>> &distribution,
      double when,
      bool previously_enabled,
      RNG& rng) {
    transitions_.Enable(tkey, distribution, previously_enabled);
  }

  virtual void Disable(const TransitionKey& tkey, double when) override {
    transitions_.Disable(tkey);
  }

  virtual void Fire(const TransitionKey& tkey, double when, RNG& rng) override {
    transitions_.Disable(tkey);
  }

 private:
  TransitionsStoreEnabled<TransitionKey,RNG> transitions_;
  std::unique_ptr<IncludeDistributions<RNG>> include_;
};


/*! For every transition, record its cumulative internal firing time.
 *  DistributionType is the unique_ptr to a distribution.
 *  QueueHandle is a handle to an entry in a mutable queue.
 */
template<typename DistributionType,typename QueueHandle>
struct TransitionTimes
{
  DistributionType dist; // nullptr for disabled transitions.
  double remaining_exponential_interval; // internal clock for transition.
  double last_modification_time; // absolute system time.
  QueueHandle queue_iter; // Points to entry in mutable priority queue.

  TransitionTimes() : dist{nullptr}, remaining_exponential_interval{0.0},
      last_modification_time(0.0), queue_iter(nullptr) {}
  TransitionTimes(DistributionType& d, double remaining, double last_mod)
      : dist(std::move(d)), remaining_exponential_interval(remaining),
      last_modification_time(last_mod), queue_iter(nullptr)
  {}

  ~TransitionTimes() {}
  // We add these because it is The Right Thing To Do
  // and because the class cannot hold a unique_ptr unless it
  // explicitly handles move operations. dist is a unique_ptr.
  TransitionTimes(TransitionTimes&& o) : dist{std::move(o.dist)},
    remaining_exponential_interval{o.remaining_exponential_interval},
    last_modification_time{o.last_modification_time},
    queue_iter{o.queue_iter} {
  }

  TransitionTimes& operator=(TransitionTimes&& o) {
    if (*this!=o) {
      dist=std::move(o);
      remaining_exponential_interval=o.remaining_exponential_interval;
      last_modification_time=o.last_modification_time;
      queue_iter=o.queue_iter;
    }
  }

  TransitionTimes(const TransitionTimes&) = delete;
  TransitionTimes& operator=(const TransitionTimes&) = delete;
};


/*! An absolute time firing time, to put into the priority_queue.
 */
template<typename TransitionKey>
struct TransitionFiringTime
{
  TransitionKey key;
  double time;
};


/*! The priority queue has to order entries in reverse time sort.
 */
template<typename TransitionKey>
struct CompareFiringTimes
{
  bool operator()(const TransitionFiringTime<TransitionKey>& a,
      const TransitionFiringTime<TransitionKey>& b) const {
    return a.time > b.time;
  }
};


/*! Anderson's method for Next Reaction of Competing Processes.
 */
template<typename TransitionKey,typename RNG>
class NonHomogeneousPoissonProcesses
: public ContinuousPropagator<TransitionKey,RNG>
{
  using FiringEntry=TransitionFiringTime<TransitionKey>;
  using Comparator=CompareFiringTimes<TransitionKey>;
  using Queue=boost::heap::fibonacci_heap<FiringEntry,
      boost::heap::compare<Comparator>>;
  using QueueHandle=typename Queue::handle_type;
  using TransitionEntry=TransitionTimes<
      std::unique_ptr<TransitionDistribution<RNG>>,QueueHandle>;
  using TransitionMap=std::map<TransitionKey,TransitionEntry>;
 public:
  NonHomogeneousPoissonProcesses()
  : include_{new BoundedHazardDistributions<RNG>{}},
    unit_exponential_{1.0}
  {}

  /*! The Include entry decides which transitions this class
   *  will take from the core matrix.
   */
  virtual bool Include(const TransitionDistribution<RNG>& distribution
      ) const override {
    return include_->Include(distribution);
  }

  /*! What transition do we think should be next?
   */
  virtual std::tuple<TransitionKey, double> Next(double now, RNG& rng
      ) const override {
    const auto& least=queue_.top();
    return std::make_tuple(least.key, least.time);
  }

  /*! Ask whether a transition is enabled.
   */
  virtual std::tuple<bool,double> Enabled(const TransitionKey& tkey
      ) const override {
    const auto& previous_distribution=distributions_.find(tkey);
    bool previously_enabled=previous_distribution!=distributions_.end();

    if (previously_enabled) {
      const TransitionEntry &entry=previous_distribution->second;
      if (entry.dist!=nullptr) {
        double enabling_time=entry.dist->EnablingTime();
        return std::make_tuple(true, enabling_time);
      } else {
        return std::make_tuple(false, 0.0);
      }
    } else {
      return std::make_tuple(false, 0.0);
    }
  }

  /*! Enable a transition. Every transition is logically already-enabled
   *  for this sampling method, so it creates on the fly transitions
   *  that were never enabled before.
   */
  virtual void Enable(const TransitionKey& tkey,
      std::unique_ptr<TransitionDistribution<RNG>>& distribution,
      double when,
      bool previously_enabled,
      RNG& rng) override {
    typename decltype(distributions_)::iterator trans_loc;
    if (!previously_enabled) {
      // Make a new timer for a new transition.
      double interval=unit_exponential_(rng);
      double firing_time=distribution->ImplicitHazardIntegral(interval, when);
      bool added=false;
      std::tie(trans_loc, added)=
          distributions_.emplace(tkey, TransitionEntry{distribution,
          interval, when});
      auto entry_handle=queue_.push({tkey, firing_time});
      trans_loc->second.queue_iter=entry_handle;
    } else {
      // Change the predicted time for a transition.
      trans_loc=distributions_.find(tkey);
      BOOST_ASSERT_MSG(trans_loc!=distributions_.end(), "A transition was "
          "enabled but cannot be found.");
      TransitionEntry& entry=trans_loc->second;
      queue_.update(entry.queue_iter,
          {tkey, distribution->ImplicitHazardIntegral(
          entry.remaining_exponential_interval, when)});
      entry.dist=std::move(distribution);
      entry.last_modification_time=when;
    }
  }

  /*! Disabling transitions means marking them as disabled and
   *  stopping their clocks. The class will remember the internal
   *  clock time, however.
   */
  virtual void Disable(const TransitionKey& tkey, double when) override {
    TransitionEntry &trans=distributions_.at(tkey);
    trans.remaining_exponential_interval-=trans.dist->HazardIntegral(
        trans.last_modification_time, when);
    trans.last_modification_time=when;
    trans.dist.reset(nullptr);
    queue_.update(trans.queue_iter,
        {tkey, (std::numeric_limits<double>::max)()});
  }


  virtual void Fire(const TransitionKey& tkey, double when, RNG& rng) override {
    BOOST_ASSERT_MSG(queue_.top().key==tkey, "Firing a transition but its "
        "key doesn't match the top key in the queue.");
    TransitionEntry &trans=distributions_.at(tkey);
    trans.remaining_exponential_interval=unit_exponential_(rng);
    double firing_time=trans.dist->ImplicitHazardIntegral(
      trans.remaining_exponential_interval, when);
    queue_.update(trans.queue_iter, {tkey, firing_time});
    trans.last_modification_time=when;
    trans.dist.reset(nullptr);
  }

 private:
  TransitionMap distributions_;
  Queue queue_;
  std::exponential_distribution<double> unit_exponential_;
  std::unique_ptr<IncludeDistributions<RNG>> include_;
};


} // smv
} // afidd


#endif // _CONTINUOUS_DYNAMICS_H_
