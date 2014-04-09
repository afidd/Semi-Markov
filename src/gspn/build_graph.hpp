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
#ifndef _BUILD_GRAPH_H_
#define _BUILD_GRAPH_H_ 1

#include <map>
#include "boost/property_map/property_map.hpp"
#include "boost/graph/copy.hpp"
#include "boost/utility/value_init.hpp"
#include "logging.hpp"
#include "petri_graph.hpp"



namespace afidd
{
namespace smv
{

/*! Maintains a mapping between PlaceKey, TransitionKey, and vertex_descriptor.
 */
template<typename BGPlace, typename BGTransition, typename vert_t>
struct BiGraphCorrespondence
{
  typedef vert_t vertex_descriptor;

  std::map<BGPlace,vert_t> pv;
  std::map<vert_t,BGPlace> vp;
  std::map<BGTransition,vert_t> tv;
  std::map<vert_t,BGTransition> vt;

  inline friend
  std::ostream& operator<<(std::ostream& os, const BiGraphCorrespondence& bgc)
  {
    return os << "BiGraphCorrespondence(" << bgc.pv.size() << ", "
      << bgc.vp.size() << ", " << bgc.tv.size() << ", " << bgc.vt.size() << ")";
  }
};



template<typename BGPlace, typename BGTransition, typename vert_t>
bool PutPlace(BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    vert_t k, BGPlace val)
{
  if (map.pv.find(val)!=map.pv.end()) {
    BOOST_LOG_TRIVIAL(error) << "Place "<<val<<" already exists.";
    return false;
  }
  map.pv.emplace(val, k);
  map.vp.emplace(k, val);
  return true;
}



template<typename BGPlace, typename BGTransition, typename vert_t>
bool PutTransition(BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    vert_t k, BGTransition val)
{
  if (map.tv.find(val)!=map.tv.end()) {
    BOOST_LOG_TRIVIAL(error) << "Transition "<<val<<" already exists.";
    return false;
  }
  map.tv.emplace(val, k);
  map.vt.emplace(k, val);
  return true;
}



template<typename BGPlace, typename BGTransition, typename vert_t>
typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::vertex_descriptor
GetPvertex(const BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    BGPlace& key)
{
  auto it=map.pv.find(key);
  if (it==map.pv.end()) {
    BOOST_LOG_TRIVIAL(error) << "Place does not exist: "<<key
      <<" map size " << map.pv.size();
    return vert_t{};
  }
  return it->second;
}




template<typename BGPlace, typename BGTransition, typename vert_t>
typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::vertex_descriptor
GetTvertex(const BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    BGTransition& key)
{
  auto it=map.tv.find(key);
  if (it==map.tv.end()) {
    BOOST_LOG_TRIVIAL(error) << "Transition does not exist: "<<key;
    return vert_t{};
  }
  return it->second;
}




template<typename BGPlace, typename BGTransition, typename vert_t>
BGPlace GetPlace(const BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
  typename BiGraphCorrespondence<BGPlace,BGTransition,
      vert_t>::vertex_descriptor& key)
{
  auto it=map.vp.find(key);
  if (it==map.vp.end()) {
    BOOST_LOG_TRIVIAL(error) << "Vertex for place does not exist: "<<key;
    return BGPlace{};
  }
  return it->second;
}




template<typename BGPlace, typename BGTransition, typename vert_t>
BGTransition GetTransition(
  const BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
  typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::vertex_descriptor& key)
{
  auto it=map.vt.find(key);
  if (it==map.vt.end()) {
    BOOST_LOG_TRIVIAL(error) << "Vertex for transition does not exist: "<<key;
    return BGTransition{};
  }
  return it->second;
}



/*! Build an ExplicitTransitions class.
 *  Template argument ET is the type of an ExplicitTransitions class.
 *  This uses a Boost::graph made with lists so that it preserves the
 *  order of input and output edges but doesn't need to know ahead of
 *  time how many vertices will be in the system. Then it copies that
 *  graph to a Boost::graph which has O(1) access to edges.
 */
template<typename ET>
class BuildGraph
{
  using BGPlace=typename ET::UserPlaceKey;
  using BGTransition=typename ET::UserTransitionKey;
  using Transition=typename ET::Transition;

  PetriBuildGraphType g_;
  using vert=boost::graph_traits<PetriBuildGraphType>::vertex_descriptor;
  using edge_t=boost::graph_traits<PetriBuildGraphType>::edge_descriptor;

  BiGraphCorrespondence<BGPlace,BGTransition,vert> bimap_;
  std::map<BGTransition,std::unique_ptr<Transition>> transitions_;

  // This is the translation map for the compiled graph,
  // not for the one we use internally.
  using BiMap=BiGraphCorrespondence<BGPlace,BGTransition,
      boost::graph_traits<PetriGraphType>::vertex_descriptor>;

public:
  using PlaceEdge=std::tuple<BGPlace,int>;

  BuildGraph() {}


  bool AddPlace(BGPlace p, int token_layer=0) {
    auto v=add_vertex({PetriGraphColor::Place, token_layer}, g_);
    return PutPlace(bimap_, v, p);
  }



  bool AddTransition(BGTransition t, std::vector<PlaceEdge> e,
    std::unique_ptr<Transition> transition) {
    auto tv=add_vertex({PetriGraphColor::Transition, 0}, g_);
    bool added=PutTransition(bimap_, tv, t);

    for (auto edge : e) {
      BGPlace p=std::get<0>(edge);
      int weight=std::get<1>(edge);
      if (!this->AddEdge(t, p, weight)) {
        added=false;
      }
    }

    transitions_.emplace(t, std::move(transition));

    return added;
  }



  bool AddEdge(BGTransition t, BGPlace p, int weight)
  {
    auto tv=GetTvertex(bimap_, t);
    assert(g_[tv].color==PetriGraphColor::Transition);
    auto pv=GetPvertex(bimap_, p);
    assert(g_[pv].color==PetriGraphColor::Place);

    if (weight<0) {
      edge_t new_edge;
      bool success;
      std::tie(new_edge, success)=boost::add_edge(pv, tv, {weight}, g_);
      if (!success) {
        BOOST_LOG_TRIVIAL(error) << "Could not add edge";
        return false;
      }
    } else {
      edge_t new_edge;
      bool success;
      std::tie(new_edge, success)=boost::add_edge(tv, pv, {weight}, g_);
      if (!success) {
        BOOST_LOG_TRIVIAL(error) << "Could not add edge";
        return false;
      }
    }
    return true;
  }



  /*! Translate the graph from a buildable type to a fast type.
   *  The buildable type is based on lists and has pointers
   *  for vertex_property. The fast time is based on vectors
   *  and uses int64_t as the index type.
   */
  ET Build()
  {
    // Check the current graph before we continue.
    size_t stoch_start=NumStochiometricCoefficients(g_);
    BOOST_LOG_TRIVIAL(debug)<< stoch_start
        << " stochiometric coefficients start";
    auto bipartite=IsBipartitePetriGraph(g_);
    BOOST_LOG_TRIVIAL(debug)<< " is bipartite "<<bipartite;
    assert(bipartite);

    ET et(num_vertices(g_));
    PetriGraphType& g=et.graph;
    using vert_n=boost::graph_traits<PetriGraphType>::vertex_descriptor;

    // In order to use the boost copy_graph algorithm, the original
    // graph has to have a vertex_index_t property, OR you have to 
    // make a map from vertex_descriptor to an index, which is what
    // we do here. The indices we assign here will _not_ be the same
    // as the indices in the resulting graph.
    std::map<vert,int64_t> vertex_index;
    using viter=boost::graph_traits<PetriBuildGraphType>::vertex_iterator;
    viter ind_begin;
    viter ind_end;
    std::tie(ind_begin, ind_end)=vertices(g_);
    for (int64_t idx=0; ind_begin!=ind_end; ++ind_begin, ++idx) {
      vertex_index.emplace(*ind_begin, idx);
    }
    using IndMap=boost::associative_property_map<std::map<vert,int64_t>>;
    IndMap vertex_index_map(vertex_index);

    // The orig_to_copy argument is a way to save a translation
    // table from the old vertex_descriptor to the new vertex_descriptor.
    // All examples use the vector map, which uses a vertex_index_map
    // to get int64_t offsets, but let's just use a simple map.
    // While we created a vertex_index above, it will not be the same
    // index as the new graph's vertex_descriptor, so we need this
    // translation.
    using TransMap=std::map<vert,vert_n>;
    TransMap translate;
    using TransPMap=boost::associative_property_map<TransMap>;
    TransPMap translate_pmap(translate);

    copy_graph(g_, g, boost::vertex_index_map(vertex_index_map).
      orig_to_copy(translate_pmap));

    typename ET::BiMap& b=et.bimap_;

    for (auto trans_iter=translate.begin();
        trans_iter!=translate.end();
        ++trans_iter)
    {
      auto old_vertex=trans_iter->first;
      auto new_vertex=trans_iter->second;
      // Is the old vertex a place or a transition?
      auto place_iter=bimap_.vp.find(old_vertex);
      auto transition_iter=bimap_.vt.find(old_vertex);
      if (place_iter!=bimap_.vp.end()) {
        if (transition_iter!=bimap_.vt.end()) {
          BOOST_LOG_TRIVIAL(error)<<"The same vertex points both to a place "
            "and a transition.";
          assert(transition_iter==bimap_.vt.end());
        }
        PutPlace(b, static_cast<int64_t>(new_vertex), place_iter->second);
      } else if (transition_iter!=bimap_.vt.end()) {
        PutTransition(b, static_cast<int64_t>(new_vertex),
            transition_iter->second);

        auto trans_obj_iter=transitions_.find(transition_iter->second);
        if (trans_obj_iter!=transitions_.end()) {
          et.transitions.emplace(new_vertex, std::move(trans_obj_iter->second));
        } else {
          BOOST_LOG_TRIVIAL(error) << "Every transition must have a transition"
            " object.";
          assert(trans_obj_iter!=transitions_.end());
        }
      } else {
        BOOST_LOG_TRIVIAL(error) << "A vertex wasn't translated during the "
          "build step.";
      }
    }

    // Check out whether it looks right.
    BOOST_LOG_TRIVIAL(debug)<<b;
    size_t stoch_cnt=NumStochiometricCoefficients(g);
    BOOST_LOG_TRIVIAL(debug)<< stoch_cnt << " stochiometric coefficients found";
    bool n_bipartite=IsBipartitePetriGraph(g);
    if (!n_bipartite)
    {
      BOOST_LOG_TRIVIAL(error)<< "The graph is not bipartite.";
      assert(n_bipartite);
    }
    return std::move(et);
  }
};


} // namespace smv
} // namespace afidd

#endif // _BUILD_GRAPH_H_
