#ifndef _GSPN_RANDOM_H_
#define _GSPN_RANDOM_H_ 1

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

template<typename IndexType, typename RNG>
struct BoostToUniform
{
  RNG& _rng;
  BoostToUniform(RNG& rng) : _rng(rng) {}
  IndexType operator()(IndexType cnt)
  {
    boost::random::uniform_int_distribution<IndexType> gen_idx(0, cnt-1);
    return gen_idx(_rng);
  }
};

} // end namespace detail
	


template<typename RNG>
double uniform(RNG& rng)
{
  return boost::random::uniform_01<double>()(rng);
}


template<typename RNG>
size_t uniform_index(RNG& rng, size_t cnt)
{
  boost::random::uniform_int_distribution<size_t> gen_idx(0, cnt-1);
  return gen_idx(rng);
}


template<typename RandomAccessIterator, typename RNG>
void random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
    RNG& rng)
{
  using Index=
    typename std::iterator_traits<RandomAccessIterator>::difference_type;
  detail::BoostToUniform<Index,RNG> shuffle_gen(rng);
  std::random_shuffle(first, last, shuffle_gen);
}

} // end namespace smv
} // end namespace afidd


#endif // _GSPN_RANDOM_H_
