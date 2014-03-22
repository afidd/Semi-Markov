#ifndef _BUILD_GRAPH_H_
#define _BUILD_GRAPH_H_ 1

#include <map>
#include "boost/property_map/property_map.hpp"
#include "boost/graph/copy.hpp"
#include "boost/utility/value_init.hpp"
#include "logging.h"
#include "petri_graph.h"



namespace afidd
{
  namespace smv
  {


template<typename BGPlace, typename BGTransition, typename vert_t>
struct BiGraphCorrespondence
{
  typedef vert_t vertex_descriptor;

  std::map<BGPlace,vert_t> pv;
  std::map<vert_t,BGPlace> vp;
  std::map<BGTransition,vert_t> tv;
  std::map<vert_t,BGTransition> vt;
};



template<typename BGPlace, typename BGTransition, typename vert_t>
void put_place(BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    vert_t k, BGPlace val)
{
  map.pv.emplace(val, k);
  map.vp.emplace(k, val);
}



template<typename BGPlace, typename BGTransition, typename vert_t>
void put_transition(BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    vert_t k, BGTransition val)
{
  map.tv.emplace(val, k);
  map.vt.emplace(k, val);
}


template<typename BGPlace, typename BGTransition, typename vert_t>
typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::vertex_descriptor
get_place(BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    BGPlace& key)
{
  auto it=map.pv.find(key);
  if (it==map.pv.end())
  {
    BOOST_LOG_TRIVIAL(error) << "Place does not exist: "<<key;
    return vert_t{};
  }
  return it->second;
}




template<typename BGPlace, typename BGTransition, typename vert_t>
typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::vertex_descriptor
get_transition(BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    BGTransition& key)
{
  auto it=map.tv.find(key);
  if (it==map.tv.end())
  {
    BOOST_LOG_TRIVIAL(error) << "Transition does not exist: "<<key;
    return vert_t{};
  }
  return it->second;
}





template<typename BGPlace, typename BGTransition>
class BuildGraph
{
  PetriBuildGraphType _g;
  using vert=boost::graph_traits<PetriBuildGraphType>::vertex_descriptor;
  using edge_t=boost::graph_traits<PetriBuildGraphType>::edge_descriptor;

  BiGraphCorrespondence<BGPlace,BGTransition,vert> _bimap;

public:
  // This is the translation map for the compiled graph,
  // not for the one we use internally.
  using BiMap=BiGraphCorrespondence<BGPlace,BGTransition,
      boost::graph_traits<PetriGraphType>::vertex_descriptor>;
  using PlaceEdge=std::tuple<BGPlace,int>;

  BuildGraph()
  {}


  bool add_place(BGPlace p, size_t token_layer)
  {
    auto v=add_vertex({PetriGraphColor::Place, token_layer}, _g);
    put_place(_bimap, v, p);
    return true;
  }



  bool add_transition(BGTransition t, std::vector<PlaceEdge> e)
  {
    auto tv=add_vertex({PetriGraphColor::Transition, 0}, _g);
    put_transition(_bimap, tv, t);

    for (auto edge : e)
    {
      BGPlace p=std::get<0>(edge);
      int weight=std::get<1>(edge);
      this->add_edge(t, p, weight);
    }
    return true;
  }



  bool add_edge(BGTransition t, BGPlace p, int weight)
  {
    auto tv=get_transition(_bimap, t);
    auto pv=get_place(_bimap, p);

    if (weight<0)
    {
      edge_t new_edge;
      bool success;
      std::tie(new_edge, success)=boost::add_edge(pv, tv, {weight}, _g);
      if (!success)
      {
        BOOST_LOG_TRIVIAL(error) << "Could not add edge";
        return false;
      }
    }
    else
    {
      edge_t new_edge;
      bool success;
      std::tie(new_edge, success)=boost::add_edge(tv, pv, {weight}, _g);
      if (!success)
      {
        BOOST_LOG_TRIVIAL(error) << "Could not add edge";
        return false;
      }
    }
    return true;
  }



  /*! Translate the graph from a buildable type to a fast type.
   *  The buildable type is based on lists and has pointers
   *  for vertex_property. The fast time is based on vectors
   *  and uses size_t as the index type.
   */
  std::tuple<PetriGraphType,BiMap> compile()
  {
    PetriGraphType g(num_vertices(_g));
    using vert_n=boost::graph_traits<PetriGraphType>::vertex_descriptor;

    // The index map turns vertex_descriptor from the build graph into
    // an integer index into the list of vertex_descriptors.
    //using index_map_t=boost::property_map<PetriBuildGraphType,
    //    boost::vertex_index_t>::type;
    //using IsoMap=boost::iterator_property_map<typename std::vector<vert>::iterator,
    //    index_map_t,vert_n,vert_n&>;

    std::map<vert,size_t> vertex_index;
    using viter=boost::graph_traits<PetriBuildGraphType>::vertex_iterator;
    viter ind_begin;
    viter ind_end;
    std::tie(ind_begin, ind_end)=vertices(_g);
    for (size_t idx=0; ind_begin!=ind_end; ++ind_begin, ++idx)
    {
      vertex_index.emplace(*ind_begin, idx);
    }
    using IndMap=boost::associative_property_map<std::map<vert,size_t>>;
    IndMap vertex_index_map(vertex_index);

    using TransMap=std::map<vert,vert_n>;
    TransMap translate;
    using TransPMap=boost::associative_property_map<TransMap>;
    TransPMap translate_pmap(translate);

    using TransVec=std::vector<vert_n>;
    TransVec trans_vec(num_vertices(_g));
    using TransPVec=boost::iterator_property_map<TransVec::iterator,
        IndMap,vert_n,vert_n&>;
    TransPVec translate_pvec(trans_vec.begin(), vertex_index_map);

    //std::vector<vert> iso_vec(num_vertices(g));
    //IsoMap iso_map(iso_vec.begin(), get(boost::vertex_index, _g));
    copy_graph(_g, g, boost::vertex_index_map(vertex_index_map).
      orig_to_copy(translate_pmap));


    BiMap b;
    for (auto trans_iter=translate.begin();
        trans_iter!=translate.end();
        ++trans_iter)
    {
      auto old_vertex=trans_iter->first;
      auto new_vertex=trans_iter->second;
      // Is the old vertex a place or a transition?
      auto place_iter=_bimap.vp.find(old_vertex);
      auto transition_iter=_bimap.vt.find(old_vertex);
      if (place_iter!=_bimap.vp.end())
      {
        put_place(b, new_vertex, place_iter->second);
      }
      else if (transition_iter!=_bimap.vt.end())
      {
        put_transition(b, new_vertex, transition_iter->second);
      }
      else
      {
        BOOST_LOG_TRIVIAL(error) << "A vertex wasn't translated";
      }
    }
    return std::make_tuple(g, b);
  }
};


} // namespace smv
} // namespace afidd

#endif // _BUILD_GRAPH_H_
