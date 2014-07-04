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
#include "gspn_random.hpp"
#include "local_marking.hpp"


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
  static constexpr int _layer_cnt=boost::mpl::size<container_types>::value;
//private:
  std::set<Place> modified_;
  Maps maps_;

public:
  typedef LocalMarking<Ts...> LocalTyped;

  Marking() {}
  const std::set<place_t>& Modified() const { return modified_; }
  void Clear() { modified_.clear(); }


  LocalTyped GetLocalMarking()
  {
    return LocalTyped{};
  }


  void InitLocal(LocalMarking<Ts...>& lm,
      std::vector<std::tuple<place_t,int,int>> neighbor_places) {
    lm.Resize(neighbor_places.size());
    int idx=0;
    detail::initialize_local<_layer_cnt,Maps,place_t,LocalTyped> il;

    for (auto& line : neighbor_places) {
      auto& place_id=std::get<0>(line);
      auto layer=std::get<1>(line);
      auto stochiometric_coefficient=std::get<2>(line);

      SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"Marking::init_local<"<<layer<<">("
        <<place_id<<", "<<idx<<", "<<stochiometric_coefficient<<")");
      il(maps_, place_id, idx, layer, stochiometric_coefficient, lm);

      ++idx;
    }
  }




  void ReadLocal(const LocalMarking<Ts...>& lm,
    std::vector<std::tuple<place_t,int,int>> neighbor_places) {
    const auto changes=lm.Changes();

    // Modified places just need to be put on the modified list.
    for (int midx : changes[0]) {
      modified_.insert(std::get<0>(neighbor_places.at(midx)));
    }

    // Added places need to be copied here.
    detail::add_by_layer<_layer_cnt,Maps,place_t,LocalTyped> abl;
    for (int aidx : changes[1]) {
      auto& neighbor_line=neighbor_places.at(aidx);
      auto& place_id=std::get<0>(neighbor_line);
      auto layer=std::get<1>(neighbor_line);

      abl(maps_, place_id, layer, lm, aidx);
      modified_.insert(place_id);
    }

    // Removed places need to be deleted.
    detail::erase_by_layer<_layer_cnt,Maps,place_t> eraser;
    for (int ridx : changes[2]) {
      auto& neighbor_line=neighbor_places.at(ridx);
      auto& place_id=std::get<0>(neighbor_line);
      auto layer=std::get<1>(neighbor_line);

      eraser(maps_, place_id, layer);
      modified_.insert(place_id);
    }
  }



  inline friend
  std::ostream& operator<<(std::ostream& os, const Marking& m) {
    const auto& mmap=std::get<0>(m.maps_);
    for (auto& kv : mmap) {
      os << "("<< kv.first<<","<<kv.second.size() <<") ";
    }
    return os;
  }
};





template<int I, typename Marking>
void Add(Marking& m, typename Marking::place_t place_id,
  const typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type& token) {
  typedef typename boost::mpl::at<typename Marking::container_types,
    boost::mpl::int_<I>>::type container_type;

  auto& typed_dict=std::get<I>(m.maps_);
  auto place_tokens=typed_dict.find(place_id);
  int added=0;
  if (place_tokens!=typed_dict.end()) {
    detail::add_to_container(place_tokens->second, token);
  } else {
    container_type c;
    detail::add_to_container(c, token);
    typed_dict.emplace(place_id, c);
  }

  m.modified_.insert(place_id);
}



template<int I, typename Marking, typename RNG>
void Remove(Marking& m, typename Marking::place_t place_id,
  size_t cnt, RNG& rng) {
  typedef typename boost::mpl::at<typename Marking::container_types,
    boost::mpl::int_<I>>::type container_type;

  auto& typed_dict=std::get<I>(m.maps_);
  auto place_tokens=typed_dict.find(place_id);
  int added=0;
  if (place_tokens!=typed_dict.end()) {
    detail::remove_from_container(place_tokens->second, cnt, rng);
    if (place_tokens->second.size()==0) {
      typed_dict.erase(place_id);
    }
  } else {
    // If we try to remove a token, there better be a token.
    assert(place_tokens!=typed_dict.end());
  }

  m.modified_.insert(place_id);
}





template<int I, typename Marking>
int Length(const Marking& m, typename Marking::place_t place_id)
{
  const auto& typed_dict=std::get<I>(m.maps_);
  auto place_tokens=typed_dict.find(place_id);
  if (place_tokens!=typed_dict.end()) {
    return place_tokens->second.size();
  } else {
    return 0;
  }
}




template<int I, typename Marking>
size_t Length(const Marking& m, typename Marking::place_t place_id,
  typename color_type<typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type>::type color) {
  const auto& typed_dict=std::get<I>(m.maps_);
  auto place_tokens=typed_dict.find(place_id);
  if (place_tokens!=typed_dict.end()) {
    const auto& colored=place_tokens->second;
    return (colored.find(color)!=colored.end()) ? 1 : 0;
  } else {
    return 0;
  }
}


/*! Runs a functor on a place marking, returning results.
 *  Given a functor which returns return_type, the signature is
 *  pair<return_type,bool>=get<0>(Marking, place_id, functor);
 *  The template argument is index of the type of token in the Marking.
 *  Runs only against the first token found.
 */
template<int I, typename Marking, typename UnaryOperator>
std::pair<
typename std::result_of<UnaryOperator(
  typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type
  )>::type,bool>
Get(const Marking& m, typename Marking::place_t place_id,
    const UnaryOperator& op) {
  typedef typename std::result_of<UnaryOperator(
  typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type
  )>::type return_type;

  const auto& typed_dict=std::get<I>(m.maps_);
  auto place_tokens=typed_dict.find(place_id);
  if (place_tokens!=typed_dict.end()) {
    auto begin=place_tokens->second.begin();
    if (begin!=place_tokens->second.end()) {
      return {op(*begin), true};
    } else {
      return {return_type(), false};
    }
  } else {
    return {return_type(), false};
  }
}


/*! Move from one container to another.
 *  Returns the count of how many added or removed.
 *  Can move from one type to another as long as tokens are the same.
 */


template<int I, int J, typename Marking, typename Modifier>
void Move(Marking& m, typename Marking::place_t place_from,
    typename Marking::place_t place_to, size_t cnt,
    const Modifier& modify_token) {
  SMVLOG(BOOST_LOG_TRIVIAL(trace)<< "Moving "<<cnt<<" tokens from "<<place_from
    <<" to "<<place_to);
  if (0==cnt) return;

  typedef typename boost::mpl::at<typename Marking::container_types,
    boost::mpl::int_<J>>::type container_type;
    
  auto& typed_dict=std::get<J>(m.maps_);
  auto place_tokens=typed_dict.find(place_from);
  auto dest_tokens=typed_dict.find(place_to);
  if (place_tokens!=typed_dict.end()) {
    if (dest_tokens==typed_dict.end()) {
      bool success;
      std::tie(dest_tokens, success)=
        typed_dict.emplace(place_to, container_type());
      assert(success);
    }
    for (auto didx=cnt; didx>0; --didx) {
      auto begin=place_tokens->second.begin();
      if (begin!=place_tokens->second.end()) {
        SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"about to apply");
        detail::apply_token_function(*begin, modify_token);
        SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"about to add");
        detail::add_to_container(dest_tokens->second, *begin);
        SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"about to erase");
        place_tokens->second.erase(begin);
        SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"erased");
      } else {
        // not enough tokens.
        assert(begin!=place_tokens->second.end());
      }
    }
  } else {
    assert(place_tokens!=typed_dict.end());
  }
  if (place_tokens->second.size()==0) {
    typed_dict.erase(place_tokens);
  }

  m.modified_.insert(place_from);
  m.modified_.insert(place_to);
  SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"move(Marking) modified "<<place_from
    <<","<<place_to);
  return;
}


template<int I, int J, typename Marking>
void Move(Marking& m, typename Marking::place_t place_from,
    typename Marking::place_t place_to, size_t cnt) {
  using TokenType=typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type&;
  detail::DoNothing<TokenType> nothing;
  afidd::smv::Move<I,J,Marking,detail::DoNothing<TokenType>>(
    m, place_from, place_to, cnt, nothing);
}


} // namespace smv
} // namespace afidd


#endif

