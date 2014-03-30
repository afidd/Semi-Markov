#ifndef _MARKING_H_
#define _MARKING_H_ 1


#include <vector>
#include <tuple>
#include <map>
#include <array>
#include <set>
#include <utility>
#include <algorithm>
#include <type_traits>
#include "boost/mpl/vector.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/at.hpp"
#include "boost/mpl/if.hpp"
#include "boost/mpl/fold.hpp"
#include "boost/mpl/range_c.hpp"
#include "boost/mpl/transform.hpp"
#include "boost/mpl/size.hpp"
#include "boost/utility/value_init.hpp"
#include "gspn_random.h"
#include "local_marking.h"


namespace afidd
{
namespace smv
{


/*! Marking contains a tuple of maps<Place,Container>
 *  where Container is either Colored or Uncolored and its
 *  value_type is a Token.
 */
template<typename Place,typename... Ts>
class Marking
{

public:
  typedef Place place_t;
  typedef typename detail::ContainerVector<Ts...>::type container_types;
  typedef typename detail::MarkTokens<Ts...>::type token_types;
  using Maps=typename detail::MapContainers<Place,Ts...>::result;
  static constexpr size_t _layer_cnt=boost::mpl::size<container_types>::value;
//private:
  std::set<Place> _modified;
  Maps _maps;

public:
  typedef LocalMarking<Ts...> LocalTyped;

  Marking() {}
  const std::set<place_t>& modified() const { return _modified; }
  void clear() { _modified.clear(); }


  LocalTyped local_marking()
  {
    return LocalTyped{};
  }


  void init_local(LocalMarking<Ts...>& lm,
      std::vector<std::tuple<size_t,size_t,int>> neighbor_places)
  {
    lm.resize(neighbor_places.size());
    size_t idx=0;
    detail::initialize_local<_layer_cnt,Maps,place_t,LocalTyped> il;

    for (auto& line : neighbor_places)
    {
      auto& place_id=std::get<0>(line);
      auto layer=std::get<1>(line);
      auto stochiometric_coefficient=std::get<2>(line);

      il(_maps, place_id, idx, layer, stochiometric_coefficient, lm);

      ++idx;
    }
  }




  void read_local(const LocalMarking<Ts...>& lm,
    std::vector<std::tuple<size_t,size_t,int>> neighbor_places)
  {
    const auto changes=lm.changes();

    // Modified places just need to be put on the modified list.
    for (size_t midx : changes[0])
    {
      _modified.insert(std::get<0>(neighbor_places.at(midx)));
    }

    // Removed places need to be deleted.
    detail::erase_by_layer<_layer_cnt,Maps,place_t> eraser;
    for (size_t ridx : changes[1])
    {
      auto& neighbor_line=neighbor_places.at(ridx);
      auto& place_id=std::get<0>(neighbor_line);
      auto layer=std::get<1>(neighbor_line);

      eraser(_maps, place_id, layer);
      _modified.insert(place_id);
    }

    // Added places need to be copied here.
    detail::add_by_layer<_layer_cnt,Maps,place_t,LocalTyped> abl;
    for (size_t aidx : changes[2])
    {
      auto& neighbor_line=neighbor_places.at(aidx);
      auto& place_id=std::get<0>(neighbor_line);
      auto layer=std::get<1>(neighbor_line);

      abl(_maps, place_id, layer, lm, aidx);
      _modified.insert(place_id);
    }
  }



  inline friend std::ostream&
  operator<<(std::ostream& os, const Marking& m)
  {
    const auto& mmap=std::get<0>(m._maps);
    for (auto& kv : mmap)
    {
      os << "("<< kv.first<<","<<kv.second.size() <<") ";
    }
    return os;
  }
};





template<size_t I, typename Marking>
void add(Marking& m, typename Marking::place_t place_id,
  const typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type& token)
{
  typedef typename boost::mpl::at<typename Marking::container_types,
    boost::mpl::int_<I>>::type container_type;

  auto& typed_dict=std::get<I>(m._maps);
  auto place_tokens=typed_dict.find(place_id);
  int added=0;
  if (place_tokens!=typed_dict.end())
  {
    detail::add_to_container(place_tokens->second, token);
  }
  else
  {
    container_type c;
    detail::add_to_container(c, token);
    typed_dict.emplace(place_id, c);
  }

  m._modified.insert(place_id);
}



template<size_t I, typename Marking, typename RNG>
void remove(Marking& m, typename Marking::place_t place_id,
  size_t cnt, RNG& rng)
{
  typedef typename boost::mpl::at<typename Marking::container_types,
    boost::mpl::int_<I>>::type container_type;

  auto& typed_dict=std::get<I>(m._maps);
  auto place_tokens=typed_dict.find(place_id);
  int added=0;
  if (place_tokens!=typed_dict.end())
  {
    detail::remove_from_container(place_tokens->second, cnt, rng);
    if (place_tokens->second.size()==0)
    {
      typed_dict.erase(place_id);
    }
  }
  else
  {
    // If we try to remove a token, there better be a token.
    assert(place_tokens!=typed_dict.end());
  }

  m._modified.insert(place_id);
}





template<size_t I, typename Marking>
size_t length(const Marking& m, typename Marking::place_t place_id)
{
  const auto& typed_dict=std::get<I>(m._maps);
  auto place_tokens=typed_dict.find(place_id);
  if (place_tokens!=typed_dict.end())
  {
    return place_tokens->second.size();
  }
  else
  {
    return 0;
  }
}




template<size_t I, typename Marking>
size_t length(const Marking& m, typename Marking::place_t place_id,
  typename color_type<typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type>::type color)
{
  const auto& typed_dict=std::get<I>(m._maps);
  auto place_tokens=typed_dict.find(place_id);
  if (place_tokens!=typed_dict.end())
  {
    const auto& colored=place_tokens->second;
    return (colored.find(color)!=colored.end()) ? 1 : 0;
  }
  else
  {
    return 0;
  }
}


/*! Runs a functor on a place marking, returning results.
 *  Given a functor which returns return_type, the signature is
 *  pair<return_type,bool>=get<0>(Marking, place_id, functor);
 *  The template argument is index of the type of token in the Marking.
 *  Runs only against the first token found.
 */
template<size_t I, typename Marking, typename UnaryOperator>
std::pair<
typename std::result_of<UnaryOperator(
  typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type
  )>::type,bool>
get(const Marking& m, typename Marking::place_t place_id,
    const UnaryOperator& op)
{
  typedef typename std::result_of<UnaryOperator(
  typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type
  )>::type return_type;

  const auto& typed_dict=std::get<I>(m._maps);
  auto place_tokens=typed_dict.find(place_id);
  if (place_tokens!=typed_dict.end())
  {
    auto begin=place_tokens->second.begin();
    if (begin!=place_tokens->second.end())
    {
      return {op(*begin), true};
    }
    else
    {
      return {return_type(), false};
    }
  }
  else
  {
    return {return_type(), false};
  }
}


/*! Move from one container to another.
 *  Returns the count of how many added or removed.
 *  Can move from one type to another as long as tokens are the same.
 */


template<size_t I, size_t J, typename Marking, typename Modifier>
void
move(Marking& m, typename Marking::place_t place_from,
    typename Marking::place_t place_to, size_t cnt,
    const Modifier& modify_token)
{
  BOOST_LOG_TRIVIAL(trace)<< "Moving "<<cnt<<" tokens from "<<place_from
    <<" to "<<place_to;
  if (0==cnt) return;

  typedef typename boost::mpl::at<typename Marking::container_types,
    boost::mpl::int_<J>>::type container_type;
    
  auto& typed_dict=std::get<J>(m._maps);
  auto place_tokens=typed_dict.find(place_from);
  auto dest_tokens=typed_dict.find(place_to);
  if (place_tokens!=typed_dict.end())
  {
    if (dest_tokens==typed_dict.end())
    {
      bool success;
      std::tie(dest_tokens, success)=
        typed_dict.emplace(place_to, container_type());
      assert(success);
    }
    for (auto didx=cnt; didx>0; --didx)
    {
      auto begin=place_tokens->second.begin();
      if (begin!=place_tokens->second.end())
      {
        BOOST_LOG_TRIVIAL(trace)<<"about to apply";
        detail::apply_token_function(*begin, modify_token);
        BOOST_LOG_TRIVIAL(trace)<<"about to add";
        detail::add_to_container(dest_tokens->second, *begin);
        BOOST_LOG_TRIVIAL(trace)<<"about to erase";
        place_tokens->second.erase(begin);
        BOOST_LOG_TRIVIAL(trace)<<"erased";
      }
      else
      {
        // not enough tokens.
        assert(begin!=place_tokens->second.end());
      }
    }
  }
  else
  {
    assert(place_tokens!=typed_dict.end());
  }
  if (place_tokens->second.size()==0)
  {
    typed_dict.erase(place_tokens);
  }

  m._modified.insert(place_from);
  m._modified.insert(place_to);
  BOOST_LOG_TRIVIAL(trace)<<"move(Marking) modified "<<place_from<<","<<place_to;
  return;
}


template<size_t I, size_t J, typename Marking>
void
move(Marking& m, typename Marking::place_t place_from,
    typename Marking::place_t place_to, size_t cnt)
{
  using TokenType=typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type&;
  detail::DoNothing<TokenType> nothing;
  afidd::smv::move<I,J,Marking,detail::DoNothing<TokenType>>(
    m, place_from, place_to, cnt, nothing);
}




} // namespace smv
} // namespace afidd


#endif

