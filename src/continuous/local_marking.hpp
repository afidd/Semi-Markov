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
#ifndef _LOCAL_MARKING_H_
#define _LOCAL_MARKING_H_ 1

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
#include "logging.hpp"



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

namespace mpl=boost::mpl;

// Given a place and a token, construct a map.
template<typename Place,typename T>
struct MakeMap
{
  typedef std::map<Place,T> type;
};


// Turn a variadic template into a tuple of values.
// Have to forward declare so that we can specialize.
template<typename T, typename R>
struct Tupleize {};


// Add a type to a tuple of types from a variadic template.
template<typename... Ts, typename R>
struct Tupleize<std::tuple<Ts...>,R>
{
  typedef std::tuple<Ts...,R> type;
};



// Turn a variadic template into a tuple of pointers to values.
// Have to forward declare so that we can specialize.
template<typename T, typename R>
struct PTupleize {};



template<typename... Ts, typename R>
struct PTupleize<std::tuple<Ts...>,R>
{
  typedef std::tuple<Ts..., typename std::add_pointer<R>::type> type;
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
    boost::mpl::transform<Containers,
                          typename MakeMap<Place,mpl_::_1>::type>::type Mapped;
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
      std::shuffle(ins.begin(), ins.end(), rng);
    }
    else
    {
      std::shuffle(outs.begin(), outs.end(), rng);
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



// Helps a full Marking find a container.
// Then it puts a reference to that container into the LocalMarking.
template<size_t I, typename Maps, typename PlaceKey, typename LM>
struct initialize_local
{
  void operator()(Maps& maps, PlaceKey& place_id, size_t idx,
    size_t layer, int stochiometric_coefficient, LM& lm)
  {
    static constexpr size_t J=I-1ul;
    if (layer==J)
    {
      auto& typed_dict=std::get<J>(maps);
      auto place_tokens=typed_dict.find(place_id);
      if (place_tokens!=typed_dict.end())
      {
        BOOST_LOG_TRIVIAL(trace)<<"initialize_local<"<<layer<<"> "<<place_id
          <<" "<<idx<<" "<<place_tokens->second.size();
        lm.template Set<J>(idx, &place_tokens->second,
                           stochiometric_coefficient);
      }
      else
      {
        BOOST_LOG_TRIVIAL(trace)<<"initialize_local<"<<layer<<"> "<<place_id
          <<" "<<idx<<" null";
        lm.template Set<J>(idx, nullptr, stochiometric_coefficient);
      }
    }

    initialize_local<J,Maps,PlaceKey,LM> il;
    il(maps, place_id, idx, layer, stochiometric_coefficient, lm);
  }
};



// Helps a full Marking find a container.
// Then it puts a reference to that container into the LocalMarking.
template<typename Maps, typename PlaceKey, typename LM>
struct initialize_local<0,Maps,PlaceKey,LM>
{
  void operator()(Maps& maps, PlaceKey& place_id, size_t idx,
    size_t layer, int stochiometric_coefficient, LM& lm)
  {}
};


template<size_t I, typename PTuple>
struct free_layer
{
  void operator()(PTuple& tuple_container, size_t layer)
  {
    if (I-1ul==layer) {
      delete std::get<I-1ul>(tuple_container);
    }
    free_layer<I-1ul,PTuple> sub;
    sub(tuple_container, layer);
  }
};



// The terminating call of free_layer.
template<typename PTuple>
struct free_layer<0,PTuple>
{
  void operator()(PTuple& tuple_container, size_t layer)
  {}
};




template<size_t I, typename Maps, typename PlaceKey>
struct erase_by_layer
{
  void operator()(Maps& maps, PlaceKey& place_id, size_t layer) {
    if (I-1ul==layer) {
      auto& typed_dict=std::get<I-1ul>(maps);
      auto place_tokens=typed_dict.find(place_id);
      if (place_tokens!=typed_dict.end()) {
        typed_dict.erase(place_tokens); 
      } else {
        // This could happen if a transition moves tokens into
        // a place and then erases in the same firing.
        BOOST_LOG_TRIVIAL(error)<<"Erasing a container because it "
          "was just emptied, but it does not exist.";
      }
    }

    erase_by_layer<I-1ul,Maps,PlaceKey> ebl;
    ebl(maps, place_id, layer);
  }
};





template<typename Maps, typename PlaceKey>
struct erase_by_layer<0,Maps,PlaceKey>
{
  void operator()(Maps& maps, PlaceKey& place_id, size_t layer)
  {}
};





template<size_t I, typename Maps, typename PlaceKey, typename LM>
struct add_by_layer
{
  void operator()(Maps& maps, PlaceKey& place_id, size_t layer,
    const LM& lm, size_t edge) {
    if (I-1ul==layer) {
      auto& typed_dict=std::get<I-1ul>(maps);
      auto place_tokens=typed_dict.find(place_id);
      if (place_tokens==typed_dict.end()) {
        using ContainerType=decltype(place_tokens->second);
        lm.template Get<I-1>(edge,
          [&typed_dict,&place_id](const ContainerType& ct)->int {
            typed_dict[place_id]=ct;
            return 0;
          });
      } else {
        BOOST_LOG_TRIVIAL(error)<<"The local marking claims this container "
          "isn't in the marking, but we found it. "<<place_id<<" "<<layer
          <<" "<<edge;
      }
    }

    add_by_layer<I-1ul,Maps,PlaceKey,LM> abl;
    abl(maps, place_id, layer, lm, edge);
  }
};





template<typename Maps, typename PlaceKey, typename LM>
struct add_by_layer<0,Maps,PlaceKey,LM>
{
  void operator()(Maps& maps, PlaceKey& place_id, size_t layer,
    const LM& lm, size_t edge)
  {}
};

} // namespace detail




/*! Input and output marking of a transition, an ordered list.
 *  The Tokens parameters are the types of tokens in the marking.
 *  This doesn't depend on the types of the places.
 *  This class exists in order to shield transitions from how
 *  the marking is stored. The transition accesses its edges
 *  by ordered index, not by the place, or a key to the place.
 */
template<typename... Tokens>
class LocalMarking
{
  // A tuple of pointers to TokenContainers.
  typedef typename detail::ContainerVector<Tokens...>::type container_types;
  typedef typename detail::MarkTokens<Tokens...>::type token_types;

  //typedef typename boost::mpl::transform<container_types,
  //  typename std::add_pointer<mpl_::_1>::type>::type WithPointer;

  typedef typename boost::mpl::fold<container_types, std::tuple<>,
      detail::PTupleize<mpl_::_1,mpl_::_2>>::type ptuple;
  // The core data structure is a vector which containes stochiometric
  // coefficent, token_layer, and a tuple of
  // references to the different kinds of containers of tokens, one
  // kind for each possible token_layer. So it's vector of
  // tuple< tuple<container0&,container1&>, layer, stoch_coeff>.
  using PlaceVec=std::vector<std::tuple<ptuple,size_t,int>>;
  PlaceVec _m;
  std::set<size_t> _modified;
  std::set<size_t> _added;
  std::set<size_t> _removed;

public:
  /*! Constructor. Argument is the number of inputs and outputs.
   */ 
  LocalMarking() : _m{}
  {}

  ~LocalMarking() {
    for (auto ridx : _added) {
      auto& place_container=_m.at(ridx);
      auto& token_container=std::get<0>(place_container);
      auto layer=std::get<1>(place_container);

      detail::free_layer<std::tuple_size<ptuple>::value,ptuple> f{};
      f(token_container, layer);
    }
  }


  void Reserve(size_t cnt)
  {
    BOOST_LOG_TRIVIAL(trace)<<"LocalMarking::reserve "<<cnt;
    _m.reserve(cnt);
  }


  void Resize(size_t cnt)
  {
    _m.resize(cnt);
  }


  std::array<std::set<size_t>,3> Changes() const
  {
    return {_modified, _added, _removed};
  }


  /*! Set pointers to the tokens at each place. Initialization of this object.
   */
  template<size_t I>
  void Set(size_t place_idx, typename boost::mpl::at<
    typename LocalMarking::container_types,boost::mpl::int_<I>>::type* ptr,
    int stochiometric_coefficient) {
    auto& place_container=_m.at(place_idx);

    auto& token_container=std::get<0>(place_container);
    std::get<I>(token_container)=ptr;

    std::get<1>(place_container)=I;
    std::get<2>(place_container)=stochiometric_coefficient;
  }




  template<size_t I>
  void Add(size_t place_idx, const typename boost::mpl::at<
      typename LocalMarking::token_types,boost::mpl::int_<I>>::type& token)
  {
    BOOST_LOG_TRIVIAL(trace)<<"LocalMarking::add";
    typedef typename boost::mpl::at<typename LocalMarking::container_types,
      boost::mpl::int_<I>>::type container_type;

    auto& place_container=_m.at(place_idx);
    auto& token_containers=std::get<0>(place_container);
    auto& token_container=std::get<I>(token_containers);

    if (token_container!=nullptr) {
      detail::add_to_container(*token_container, token);
    } else {
      token_container=new container_type{};
      _added.insert(place_idx);
      detail::add_to_container(*token_container, token);
    }

    _modified.insert(place_idx);
  }




  template<size_t I, typename RNG>
  void Remove(size_t place_idx, size_t cnt, RNG& rng) {
    BOOST_LOG_TRIVIAL(trace)<<"LocalMarking::remove";
    typedef typename boost::mpl::at<typename LocalMarking::container_types,
      boost::mpl::int_<I>>::type container_type;

    auto& place_container=_m.at(place_idx);
    auto& token_containers=std::get<0>(place_container);
    auto& token_container=std::get<I>(token_containers);
    if (token_container!=nullptr) {
      detail::remove_from_container(*token_container, cnt, rng);
    } else {
      BOOST_LOG_TRIVIAL(error)<<"LocalMarking::remove Cannot remove from "
        "nonexistent container "<<place_idx<<" "<<I<<" "<<cnt;
    }
    if (token_container->size()==0) {
      _removed.insert(place_idx);
    } else {
      _modified.insert(place_idx);
    }
  }




  template<size_t I>
  size_t Length(size_t place_idx) const {
    auto& place_container=_m.at(place_idx);
    auto& token_container=std::get<I>(std::get<0>(place_container));
    if (token_container!=nullptr) {
      return token_container->size();
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
  template<size_t I, typename UnaryOperator>
  std::pair<
  typename std::result_of<UnaryOperator(
    typename boost::mpl::at<
        typename LocalMarking::container_types,boost::mpl::int_<I>>::type
    )>::type,bool>
  Get(size_t place_idx, const UnaryOperator& op) const {
    typedef typename std::result_of<UnaryOperator(
    typename boost::mpl::at<
        typename LocalMarking::container_types,boost::mpl::int_<I>>::type
    )>::type return_type;

    auto& place_container=_m.at(place_idx);
    auto& token_containers=std::get<0>(place_container);
    auto& token_container=std::get<I>(token_containers);

    if (token_container!=nullptr) {
      return {op(*token_container), true};
    } else {
      return {return_type{}, false};
    }
  }



  /*! Runs a functor on a token, returning results.
   *  Given a functor which returns return_type, the signature is
   *  pair<return_type,bool>=get<0>(Marking, place_id, functor);
   *  The template argument is index of the type of token in the Marking.
   *  Runs only against the first token found.
   */
  template<size_t I, typename UnaryOperator>
  std::pair<
  typename std::result_of<UnaryOperator(
    typename boost::mpl::at<
        typename LocalMarking::token_types,boost::mpl::int_<I>>::type
    )>::type,bool>
  GetToken(size_t place_idx, const UnaryOperator& op) const {
    typedef typename std::result_of<UnaryOperator(
    typename boost::mpl::at<
        typename LocalMarking::token_types,boost::mpl::int_<I>>::type
    )>::type return_type;

    auto& place_container=_m.at(place_idx);
    auto& token_containers=std::get<0>(place_container);
    auto& token_container=std::get<I>(token_containers);

    if (token_container!=nullptr) {
      auto begin=token_container->begin();
      if (begin!=token_container->end()) {
        return {op(*begin), true};
      } else {
        return {return_type{}, false};
      }
    } else {
      return {return_type{}, false};
    }
  }



  template<size_t I, size_t J, typename Modifier>
  void
  Move(size_t place_from, size_t place_to, size_t cnt,
      const Modifier& modify_token)
  {
    BOOST_LOG_TRIVIAL(trace)<< "Moving "<<cnt<<" tokens from "<<place_from
      <<" to "<<place_to;
    if (0==cnt) return;

    auto& place_from_container=_m.at(place_from);
    auto& from_container=std::get<I>(std::get<0>(place_from_container));

    if (from_container!=nullptr) {
      if (from_container->size()>=cnt) {
        typedef typename boost::mpl::at<typename LocalMarking::container_types,
          boost::mpl::int_<J>>::type to_container_type;
        auto& place_to_container=_m.at(place_to);
        auto& to_container=std::get<J>(std::get<0>(place_to_container));
        if (to_container==nullptr) {
          BOOST_LOG_TRIVIAL(trace)<< "Moving to_container null";
          to_container=new to_container_type{};
          _added.insert(place_to);
        } else {
          BOOST_LOG_TRIVIAL(trace)<< "Moving to_container exists";
        }
        if (from_container!=to_container) {
          for (auto didx=cnt; didx>0; --didx) {
            auto begin=from_container->begin();
            if (begin!=from_container->end()) {
              detail::apply_token_function(*begin, modify_token);
              detail::add_to_container(*to_container, *begin);
              from_container->erase(begin);
            }
          }
        } else {
          // If the two containers are different, still apply the functor.
          auto begin=from_container->begin();
          for (auto didx=cnt; didx>0; --didx) {
            detail::apply_token_function(*begin, modify_token);
            ++begin;
          }
        }
        if (from_container->size()==0) {
          BOOST_LOG_TRIVIAL(trace)<<"Moving mark "<<place_from<<" to erase.";
          _removed.insert(place_from);
        }
      } else {
        BOOST_LOG_TRIVIAL(error)<<"Not enough tokens to move "
          <<place_from<<" "<<place_to<<" "<<I<<" "<<cnt;
      }
    } else {
      BOOST_LOG_TRIVIAL(error)<<"Cannot move a token from an empty container "
        <<place_from<<" "<<place_to<<" "<<I<<" "<<cnt;
    }
    BOOST_LOG_TRIVIAL(trace)<< "~Moving "<<cnt<<" tokens from "<<place_from;
  }



  size_t Layer(size_t place_idx) const {
    auto& place_container=_m.at(place_idx);
    return std::get<1>(place_container);
  }



  int stochiometric_coefficient(size_t place_idx) const {
    auto& place_container=_m.at(place_idx);
    return std::get<2>(place_container);
  }



  template<size_t I, size_t J>
  void Move(size_t place_from, size_t place_to, size_t cnt) {
    using TokenType=typename boost::mpl::at<
        typename LocalMarking::token_types,boost::mpl::int_<I>>::type&;

    detail::DoNothing<TokenType> nothing;

    this->template Move<I,J,detail::DoNothing<TokenType>>(
      place_from, place_to, cnt, nothing);
  }



  template<size_t I, typename RNG, typename AndModify>
  void TransferByStochiometricCoefficient(RNG& rng, const AndModify& mod) {
    BOOST_LOG_TRIVIAL(trace)<<"TransferByStochiometricCoefficient";
    using TokenType=typename boost::mpl::at<
      typename LocalMarking::token_types,boost::mpl::int_<I>>::type&;
    detail::DoNothing<TokenType> do_nothing;

    std::vector<size_t> ins;
    std::vector<size_t> outs;

    size_t place_idx=0;
    size_t layer=0;
    int weight=0;

    for (auto& collect_place : _m) {
      layer=std::get<1>(collect_place);
      weight=std::get<2>(collect_place);

      if (layer==I) {
        if (weight<0) {
          std::fill_n(std::back_inserter(ins), -weight, place_idx);
        } else if (weight>0) {
          std::fill_n(std::back_inserter(outs), weight, place_idx);
        } else {
          ; // Don't worry about stochiometric coefficients of 0.
        }
      }
      ++place_idx;
    }

    // Optional shuffling of arrays.
    if (ins.size()>outs.size()) {
      std::shuffle(ins.begin(), ins.end(), rng);
    } else {
      std::shuffle(outs.begin(), outs.end(), rng);
    }

    auto initer=ins.begin();
    auto outiter=outs.begin();
    for ( ; initer!=ins.end() && outiter!=outs.end(); ++initer, ++outiter) {
      this->template Move<I,I>(*initer, *outiter, 1ul, mod);
    }

    // If out needs extra tokens, create them.
    for ( ; outiter!=outs.end(); ++outiter) {
      using TokenType=
        typename boost::mpl::at<
          typename LocalMarking::token_types,
          boost::mpl::size_t<I>
        >::type;
      boost::value_initialized<TokenType> t;
      this->template Add<I>(*outiter, t.data());
    }

     // If in has extra tokens, destroy them.
    for ( ; initer!=ins.end(); ++initer) {
      this->Remove<I,RNG>(*initer, 1, rng);
    }
    BOOST_LOG_TRIVIAL(trace)<<"~TransferByStochiometricCoefficient";
  }



  template<size_t I, typename RNG>
  void TransferByStochiometricCoefficient(RNG& rng) {
    using TokenType=typename boost::mpl::at<
      typename LocalMarking::token_types,boost::mpl::int_<I>>::type&;
    detail::DoNothing<TokenType> do_nothing;
    
    TransferByStochiometricCoefficient<I,RNG,detail::DoNothing<TokenType>>(
      rng, do_nothing);
  }



  template<size_t I>
  bool InputTokensSufficient() const {
    BOOST_LOG_TRIVIAL(trace)<<"input_tokens_sufficient";

    size_t place_idx=0;
    size_t layer=0;
    int weight=0;

    for (auto& collect_place : _m) {
      layer=std::get<1>(collect_place);
      weight=std::get<2>(collect_place);
      if (layer==I) {
        if (weight<0) {
          auto available=static_cast<int>(this->template Length<I>(place_idx));
          if (available+weight<0) {
            return false;
          }
        } else {
          ; // Don't worry about other stochiometric coefficients.
        }
      }
      ++place_idx;
    }
    return true;
  }



  template<size_t I>
  bool OutputsTokensEmpty() const {
    size_t place_idx;
    size_t layer;
    int weight;

    for (auto collect_place : _m.place_indexes()) {
      std::tie(place_idx, layer, weight)=collect_place;
      if (layer==I) {
        if (weight>0) {
          auto available=this->template Length<I>(place_idx);
          if (available>0) {
            return false;
          }
        } else {
          ; // Don't worry about other stochiometric coefficients.
        }
      }
    }
    return true;
  }


  inline friend
  std::ostream& operator<<(std::ostream& os, const LocalMarking& lm) {
    return os << "Local marking "<<lm._m.size();
  }

};






}
}

#endif // _LOCAL_MARKING_H_
