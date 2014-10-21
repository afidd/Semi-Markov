#ifndef _ENSEMBLE_H_
#define _ENSEMBLE_H_ 1

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>


namespace afidd
{
namespace smv
{
/*! Runs a function many times over a set of random number generators.
 *  Runnable is a function that takes a random number generator and seed
 *  as argument.
 */
template<typename Runnable,typename RandGen>
class Ensemble {
  Runnable runner_;
  std::vector<RandGen> rng_;
  int thread_cnt_;
  int run_cnt_;
  size_t rand_seed_;
  // States are 0 available 1 running 2 please join this thread.
  std::vector<int> ready_flag_;
  std::vector<std::thread> thread_;
  std::mutex ensemble_m_;
  std::condition_variable thread_done_;
 public:
  Ensemble(Runnable runner, int thread_cnt, int run_cnt, size_t rand_seed)
  : runner_(runner), thread_cnt_(thread_cnt), run_cnt_(run_cnt),
  rand_seed_(rand_seed), rng_(thread_cnt), thread_(thread_cnt),
  ready_flag_(thread_cnt) {
    BOOST_LOG_TRIVIAL(info)<<"threads "<<thread_cnt<<" runs "<<run_cnt;
    for (auto& rnginit : rng_) {
      rnginit.seed(rand_seed);
      ++rand_seed;
    }
    for (size_t flag_init_idx=0; flag_init_idx<thread_cnt; ++flag_init_idx) {
      ready_flag_[flag_init_idx]=0;
    }
    BOOST_LOG_TRIVIAL(info)<<"Next available rand seed: "<<rand_seed;
  }

  void Run() {
    // Three phases: spin up, work, spin down.
    int up_thread=thread_cnt_;
    while (run_cnt_>0 && up_thread>0) {
      --up_thread;
      int seed=rand_seed_+up_thread;
      int run=run_cnt_;
      int tidx=up_thread;
      auto run_notify=[&,seed, run, tidx]()->void {
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread start "<<run<<" "<<tidx);
        runner_(rng_[tidx], seed, run);
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread finish "<<run
          <<" tidx "<<tidx);
        std::unique_lock<std::mutex> register_done(ensemble_m_);
        assert(ready_flag_[tidx]==1);
        ready_flag_[tidx]=2;
        register_done.unlock();
        thread_done_.notify_one();
      };
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread run "<<run_cnt_
        <<" up_thread "<<up_thread);
      assert(ready_flag_[up_thread]==0);
      ready_flag_[up_thread]=1;
      thread_[up_thread]=std::thread(run_notify);
      --run_cnt_;
    }

    std::unique_lock<std::mutex> checking_available(ensemble_m_);
    while (run_cnt_>0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread wait m "<<run_cnt_);
      thread_done_.wait(checking_available);
      int available=ThreadFinished();
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread avail "<<available);
      while (available>=0) {
        thread_[available].join();
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread joined "<<available);
        int seed=rand_seed_+available;
        int run=run_cnt_;
        auto run_notify=[&, seed, run, available]()->void {
          SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread start m "<<run
            <<" "<<available;);
          runner_(rng_[available], seed, run);
          SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread finish m "<<run<<
            " avail "<<available);
          std::unique_lock<std::mutex> register_done(ensemble_m_);
          ready_flag_[available]=2;
          register_done.unlock();
          thread_done_.notify_one();
        };
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread run m "<<run_cnt_
          << " avail "<<available);
        assert(ready_flag_[available]==2);
        ready_flag_[available]=1;
        thread_[available]=std::thread(run_notify);
        --run_cnt_;
        available=ThreadFinished();
      }
    }

    while (ThreadRunning()) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread wait d "<<run_cnt_);
      thread_done_.wait(checking_available);
      int available=ThreadFinished();
      while (available>=0) {
        thread_[available].join();
        assert(ready_flag_[available]==2);
        ready_flag_[available]=0;
        available=ThreadFinished();
      }
    }
  }

 private:
  int ThreadFinished() {
    int available=-1;
    for (int j=0; j<thread_cnt_; ++j) {
      if (ready_flag_[j]==2) {
        available=j;
        break;
      }
    }
    return available;
  }

  bool ThreadRunning() {
    for (auto& check_ready : ready_flag_) {
      if (check_ready!=0) return true;
    }
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"no threads running");
    return false;
  }
};

}
}
#endif
