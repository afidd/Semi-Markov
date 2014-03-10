#ifndef _MARKING_H_
#define _MARKING_H_ 1


#include <vector>
#include <tuple>
#include <map>
#include <array>
#include <set>
#include <utility>
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


namespace afidd
{
namespace smv
{
template<typename Place,typename... Ts>
class Marking;


template<typename Token>
using Uncolored=std::vector<Token>;


template<typename Token>
struct color_type
{
  typedef size_t type;
};


template<typename Token>
struct unique_color
{
  static const bool value=false;
};


template<typename Token>
const typename color_type<Token>::type color(const Token&) {
  return typename color_type<Token>::type();
};


template<typename Token>
using Colored=std::map<typename color_type<Token>::type,Token>;


namespace detail
{
using namespace std;
namespace mpl=boost::mpl;

template<typename Place,typename T>
struct MakeMap
{
  typedef std::map<Place,T> type;
};

template<typename T, typename R>
struct Tupleize {};


template<typename... Ts, typename R>
struct Tupleize<std::tuple<Ts...>,R>
{
  typedef std::tuple<Ts...,R> type;
};


// Given a vector or a map, find the type of the values
// stored, so throw out the key of the map.
template<typename Container>
struct TokenOfContainer
{
  typedef typename Container::value_type type;
};

template<typename Color, typename Token>
struct TokenOfContainer<std::map<Color,Token>>
{
  typedef Token type;
};


// A container is a Colored<Cow> or Uncolored<Cow>.
template<typename... Ts>
struct ContainerVector
{
  typedef boost::mpl::vector<Ts...> type;
};


// The map containers is std::tuple<std::map<place id, Colored<Cow>>, ...>
template<typename Place,typename... Ts>
struct MapContainers
{
  typedef boost::mpl::vector<Ts...> Containers;
  typedef typename
    boost::mpl::transform<Containers,typename MakeMap<Place,mpl_::_1>::type>::type Mapped;
  typedef typename boost::mpl::fold<Mapped, std::tuple<>,
    Tupleize<mpl_::_1,mpl_::_2>>::type result;
};


// Just the tokens, as in Cow() or Pig().
template<typename... Ts>
struct MarkTokens
{
  typedef boost::mpl::vector<Ts...> Containers;
  typedef typename
  boost::mpl::transform<Containers,TokenOfContainer<boost::mpl::_1>>::type type;
};


// This set of routines makes similar access to vectors and maps of tokens.
template<typename Container, typename Token>
int add_to_container(Container& c, const Token& t) { return 0; }


template<typename Key,typename Token>
int add_to_container(std::map<Key,Token>& c, const Token& t)
{
  c.emplace(color(t), t);
  return 1;
};


template<typename Token>
int add_to_container(std::vector<Token>& c, const Token& t)
{
  c.emplace_back(t);
  return 1;
};


template<typename Container, typename RNG>
void remove_from_container(Container& c, size_t cnt, RNG& rng)
{
  auto len=std::distance(c.begin(), c.end());
  while (cnt>0)
  {
    size_t idx=uniform_index(rng, len);
    auto to_erase=c.begin();
    std::advance(to_erase, idx);
    c.erase(to_erase);
    --cnt;
    --len;
  }
}


template<typename Token, typename Functor>
void apply_token_function(Token& token, const Functor& f)
{
  f(token);
}


template<typename Color, typename Token, typename Functor>
void apply_token_function(std::pair<Color,Token>& val, const Functor& f)
{
  f(val.second);
}


// Convenience methods to move tokens.

template<typename LM,typename RNG>
struct NoCopy
{
  static void copy_layer(LM& lm, RNG& rng) {};
};


template<typename... Args>
struct DoNothing
{
  void operator()(Args&&... args) const {}
};


template<typename First, typename Last, typename LM, typename RNG>
struct CopyTokens
{
  static void copy_layer(LM& local_mark, RNG& rng)
  {
    using LayerConstant=typename boost::mpl::deref<First>::type;
    static const size_t I=LayerConstant::value;

    std::vector<size_t> ins;
    std::vector<size_t> outs;

    typename LM::place_t place;
    size_t idx=0;
    size_t layer;
    int weight;

    for (auto collect_place : local_mark.place_indexes())
    {
      std::tie(place, layer, weight)=collect_place;
      if (layer==I)
      {
        if (weight<0)
        {
          std::fill_n(std::back_inserter(ins), -weight, idx);
        }
        else if (weight>0)
        {
          std::fill_n(std::back_inserter(outs), weight, idx);
        }
        else
        {
          ; // Don't worry about stochiometric coefficients of 0.
        }
      }
      ++idx;
    }

    // Optional shuffling of arrays.
    if (ins.size()>outs.size())
    {
      afidd::smv::random_shuffle(ins.begin(), ins.end(), rng);
    }
    else
    {
      afidd::smv::random_shuffle(outs.begin(), outs.end(), rng);
    }


    auto initer=ins.begin();
    auto outiter=outs.begin();
    for ( ; initer!=ins.end() && outiter!=outs.end(); ++initer, ++outiter)
    {
      static_assert(std::is_same<size_t,
        typename LayerConstant::value_type>::value,
        "I is not a size_t");
      local_mark.template move<I,I>(*initer, *outiter, 1);
    }

    // If out needs extra tokens, create them.
    for ( ; outiter!=outs.end(); ++outiter)
    {
      using TokenType=typename mpl::at<typename LM::Marking::token_types,LayerConstant>::type;
      boost::value_initialized<TokenType> t;
      local_mark.template add<I>(*outiter, t.data());
    }

     // If in has extra tokens, destroy them.
    for ( ; initer!=ins.end(); ++initer)
    {
      local_mark.template remove<I,RNG>(*initer, 1, rng);
    }

    typedef typename std::conditional
      <
        std::is_same<typename boost::mpl::next<First>::type, Last>::value,
        NoCopy<LM,RNG>,
        CopyTokens<typename boost::mpl::next<First>::type, Last, LM, RNG>
      >::type SubCopy;
    SubCopy::copy_layer(local_mark, rng);
  }
};



} // namespace detail



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
//private:
  std::set<Place> _modified;
  Maps _maps;
public:
  Marking() {}
  const std::set<place_t>& modified() { return _modified; }
  void clear() { _modified.clear(); }
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
  if (place_from==place_to) return;

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


////////////////////////////////////////
// LocalMarking
////////////////////////////////////////


/*! This presents that part of the marking upon which a transition depends.
 *  The places are ordered and indexed consecutively from 0.
 */
template<typename MarkingType>
class LocalMarking
{
public:
  typedef MarkingType Marking;
  typedef typename Marking::place_t place_t;
  typedef typename boost::mpl::range_c<size_t,0,boost::mpl::size<
      typename Marking::token_types>::value> layers;
  typedef LocalMarking<Marking> type;
private:
  Marking& _m;
  using PlaceVec=std::vector<std::tuple<place_t,size_t,int>>;
  // place token, token_layer, stochiometric_coefficient
  PlaceVec _places;


public:
  LocalMarking(Marking& marking, const PlaceVec& places)
  : _m(marking), _places(places) {}


  friend std::ostream& operator<<(std::ostream& os,
      const LocalMarking<MarkingType>& lm)
  {
    for (auto& t : lm._places)
    {
      os << "(" << std::get<0>(t) << "," << std::get<1>(t) << ","
        << std::get<2>(t) << ") ";
    }
    return os;
  }


  template<size_t I>
  size_t length(size_t place_idx) const
  {
    return afidd::smv::length<I>(_m, std::get<0>(_places.at(place_idx)));
  }


  const PlaceVec& place_indexes() const
  {
    return _places;
  }


  size_t layer(size_t place_idx) const
  {
    return std::get<1>(_places.at(place_idx));
  }



  int stochiometric_coefficient(size_t place_idx) const
  {
    return std::get<2>(_places.at(place_idx));
  }


/*  template<size_t I, typename... Args>
  decltype(afidd::smv::add(_m, std::forward<Args>(args)...))
  add(Args&&... args)
  -> decltype(f(std::forward<Args>(args)...))
  {
    afidd::smv::add(_m, std::forward<Args>(args)...);
  }
*/

  template<size_t I, typename... Args>
  void
  add(Args&&... args)
  {
    afidd::smv::add<I>(_m, std::forward<Args>(args)...);
  }


  template<size_t I, typename RNG>
  void remove(size_t place_idx, size_t cnt, RNG& rng)
  {
    afidd::smv::template remove<I>(_m, std::get<0>(_places.at(place_idx)), cnt, rng);
  }



  template<size_t I>
  size_t length(size_t place_idx,
    typename color_type<typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type>::type color) const
  {
    return afidd::smv::length<I>(_m, std::get<0>(_places.at(place_idx)), color);
  }



  template<size_t I, typename UnaryOperator>
  std::pair
    <
      typename std::result_of<UnaryOperator(
        typename boost::mpl::at<
            typename Marking::token_types,boost::mpl::int_<I>>::type
        )>::type,
      bool
    >
  get(size_t place_idx, const UnaryOperator& op)
  {
    auto pid=std::get<0>(_places.at(place_idx));
    return afidd::smv::get<I,Marking,UnaryOperator>(_m, pid, op);
  }



  template<size_t I, size_t J>
  void
  move(size_t place_from, size_t place_to, size_t cnt)
  {
    afidd::smv::move<I,J,Marking>(_m,
        std::get<0>(_places.at(place_from)),
        std::get<0>(_places.at(place_to)), cnt);
  }


  template<typename RNG>
  void fire_by_stochiometric_coefficient(RNG& rng)
  {
    using TopCopy=detail::CopyTokens
      <
      typename boost::mpl::begin<layers>::type,
      typename boost::mpl::end<layers>::type,
      type,
      RNG
      >;
    TopCopy::copy_layer(*this, rng);
  }



  template<size_t I, typename RNG, typename AndModify>
  void transfer_by_stochiometric_coefficient(RNG& rng, const AndModify mod)
  {
    using TokenType=typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type&;
    detail::DoNothing<TokenType> do_nothing;


    std::vector<size_t> ins;
    std::vector<size_t> outs;

    place_t place;
    size_t layer;
    int weight;

    for (auto collect_place : this->place_indexes())
    {
      std::tie(place, layer, weight)=collect_place;
      if (layer==I)
      {
        if (weight<0)
        {
          std::fill_n(std::back_inserter(ins), -weight, place);
        }
        else if (weight>0)
        {
          std::fill_n(std::back_inserter(outs), weight, place);
        }
        else
        {
          ; // Don't worry about stochiometric coefficients of 0.
        }
      }
    }

    // Optional shuffling of arrays.
    if (ins.size()>outs.size())
    {
      afidd::smv::random_shuffle(ins.begin(), ins.end(), rng);
    }
    else
    {
      afidd::smv::random_shuffle(outs.begin(), outs.end(), rng);
    }


    auto initer=ins.begin();
    auto outiter=outs.begin();
    for ( ; initer!=ins.end() && outiter!=outs.end(); ++initer, ++outiter)
    {
      afidd::smv::move<I,I>(_m, *initer, *outiter, 1, mod);
    }

    // If out needs extra tokens, create them.
    for ( ; outiter!=outs.end(); ++outiter)
    {
      using TokenType=
        typename boost::mpl::at<
          typename Marking::token_types,
          boost::mpl::size_t<I>
        >::type;
      boost::value_initialized<TokenType> t;
      this->add<I>(*outiter, t.data());
    }

     // If in has extra tokens, destroy them.
    for ( ; initer!=ins.end(); ++initer)
    {
      this->remove<I,RNG>(*initer, 1, rng);
    }

  }



  template<size_t I, typename RNG>
  void transfer_by_stochiometric_coefficient(RNG& rng)
  {
    using TokenType=typename boost::mpl::at<
      typename Marking::token_types,boost::mpl::int_<I>>::type&;
    detail::DoNothing<TokenType> do_nothing;
    
    transfer_by_stochiometric_coefficient<I,RNG,detail::DoNothing<TokenType>>(
      rng, do_nothing);
  }



  template<size_t I>
  bool input_tokens_sufficient() const
  {
    place_t place;
    size_t layer;
    int weight;

    for (auto collect_place : this->place_indexes())
    {
      std::tie(place, layer, weight)=collect_place;
      if (layer==I)
      {
        if (weight<0)
        {
          auto available=static_cast<int>(afidd::smv::length<I>(_m, place));
          if (available+weight<0)
          {
            return false;
          }
        }
        else
        {
          ; // Don't worry about other stochiometric coefficients.
        }
      }
    }
    return true;
  }



  template<size_t I>
  bool outputs_tokens_empty() const
  {
    place_t place;
    size_t layer;
    int weight;

    for (auto collect_place : _m.place_indexes())
    {
      std::tie(place, layer, weight)=collect_place;
      if (layer==I)
      {
        if (weight>0)
        {
          auto available=afidd::smv::length<I>(_m, place);
          if (available>0)
          {
            return false;
          }
        }
        else
        {
          ; // Don't worry about other stochiometric coefficients.
        }
      }
    }
    return true;
  }
};



} // namespace smv
} // namespace afidd


#endif

