#ifndef _GSPN_RANDOM_H_
#define _GSPN_RANDOM_H_ 1

#include <random>
#include <mutex>
#include "boost/random/uniform_01.hpp"
#include "boost/random/uniform_int_distribution.hpp"

/*! These routines wrap uniform generators for our code.
 *  Options include std::random or boost::random.
 *  We may want to use a lagged Fibonacci for parallel work.
 */


namespace afidd
{
namespace smv
{


namespace detail
{
/*
template<typename IndexType, typename RNG>
struct RngToUniform
{
  RNG& _rng;
  RngToUniform(RNG& rng) : _rng(rng) {}
  IndexType operator()(IndexType cnt)
  {
    return _rng()%cnt;
  }
};


template<typename IndexType, typename RNG>
struct RngToUniform<IndexType,BoostGenerator<RNG>>
{
  RNG& _rng;
  RngToUniform(BoostGenerator<RNG>& rng) : _rng(rng.rng) {}
  IndexType operator()(IndexType cnt)
  {
    boost::random::uniform_int_distribution<IndexType> gen_idx(0, cnt-1);
    return gen_idx(_rng);
  }
};



template<typename IndexType, typename RNG>
struct RngToUniform<IndexType,StandardGenerator<RNG>>
{
  RNG& _rng;
  RngToUniform(StandardGenerator<RNG>& rng) : _rng(rng.rng) {}
  IndexType operator()(IndexType cnt)
  {
    std::uniform_int_distribution<IndexType> gen_idx(0, cnt-1);
    return gen_idx(_rng);
  }
};
*/
} // end namespace detail
	

/*! Random number generator gets a list of numbers from another generator.
 *  This is a safe, simple way to get parallel random numbers. Use
 *  a single generator to make them.
 *  The template argument is a ProviderGenerator.
 */
template<typename ProviderGenerator>
class SlurpGenerator
{
 public:
  typedef unsigned long result_type;

 private:
  ProviderGenerator& _provider;
  std::vector<result_type> _cache;
  typename std::vector<result_type>::const_iterator _cur;
  size_t _idx;

 public:
  SlurpGenerator(ProviderGenerator& gen, size_t idx, size_t capacity=1000)
  : _provider(gen), _idx{idx}, _cache(capacity), _cur(_cache.begin()) {}

  result_type min() { return _provider.min(); }
  result_type max() { return _provider.max(); }



  void release()
  {
    _provider.release(_idx);
  }



  result_type operator()()
  {
    if (_cur==_cache.end())
    {
       _provider.write(_cache);
       _cur=_cache.begin();
    }
    auto res=*_cur;
    ++_cur;
    return res;
  }

};





/*! A random number generator that can spawn child generators.
 *  Only one generator is used to make all of the random numbers.
 *  The child sub-generators fill vectors of rands from the
 *  parent.
 *
 *  The child generators check out a child because current
 *  parallel codes (TBB and C++11) don't give you a thread id
 *  to index the generator. This puts some burden on the user
 *  to release generators.
 */
template<typename RandGen>
class ProviderGenerator
{
public:
  typedef typename RandGen::result_type result_type;
  typedef SlurpGenerator<ProviderGenerator> value_type;
private:
  RandGen& _gen;
  std::mutex _one_reader;
  std::vector<value_type> _subgen;
  std::set<size_t> _available;
  size_t _sub_capacity;
  size_t _too_many;

 public:
 ProviderGenerator(RandGen& gen, size_t capacity=1000,
    size_t too_many_generators=100)
  : _gen(gen), _sub_capacity(capacity), _too_many(too_many_generators)
 {
  _subgen.emplace_back(value_type(*this, 0, _sub_capacity));
 }

  result_type min() { return _gen.min(); }
  result_type max() { return _gen.max(); }

  result_type operator()()
  {
    return _subgen.at(0)();
  }


  value_type& child()
  {
    std::lock_guard<std::mutex> guard{_one_reader};
    if (_available.size()>0)
    {
      auto give_iter=_available.begin();
      auto& give_rng=_subgen.at(*give_iter);
      _available.erase(give_iter);
      return give_rng;
    }
    else
    {
      size_t idx=_subgen.size();
      if (idx>_too_many)
      {
        // Indicates the code is failing to release a generator.
        assert(idx<=_too_many);
      }
      _subgen.emplace_back(value_type(*this, idx, _sub_capacity));
      return _subgen.at(_subgen.size()-1);
    }
  }


  void release(size_t which_subgen)
  {
    std::lock_guard<std::mutex> guard{_one_reader};
    _available.insert(which_subgen);
  }


  void write(std::vector<result_type>& cache)
  {
    std::lock_guard<std::mutex> guard{_one_reader};
    for (auto& val : cache)
    {
      val=_gen();
    }
  }
};



template<typename RNG>
double uniform(RNG& rng)
{
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(rng);
}



template<typename RNG>
size_t uniform_index(RNG& rng, size_t cnt)
{
  std::uniform_int_distribution<size_t> gen_idx(0, cnt-1);
  return gen_idx(rng);
}


/*

template<typename RandomAccessIterator, typename RNG>
void random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
    RNG& rng)
{
  using Index=
    typename std::iterator_traits<RandomAccessIterator>::difference_type;
  detail::RngToUniform<Index,RNG> shuffle_gen(rng);
  std::random_shuffle(first, last, shuffle_gen);
} */

} // end namespace smv
} // end namespace afidd


#endif // _GSPN_RANDOM_H_
