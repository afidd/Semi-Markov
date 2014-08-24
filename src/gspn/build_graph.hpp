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
#include "boost/bimap.hpp"
#include "boost/bimap/unordered_set_of.hpp"
#include "boost/bimap/multiset_of.hpp"
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

  typedef boost::bimap<boost::bimaps::set_of<vert_t>,
            boost::bimaps::set_of<BGPlace>> PlaceMap;
  // There may be more than one transition with the same
  // key because we can use transition fusion to fuse nets
  // and not rename the transitions while renaming all places.
  typedef boost::bimap<boost::bimaps::set_of<vert_t>,
            boost::bimaps::multiset_of<BGTransition>> TransitionMap;

  PlaceMap pm_;
  TransitionMap tm_;

  inline friend
  std::ostream& operator<<(std::ostream& os, const BiGraphCorrespondence& bgc)
  {
    return os << "BiGraphCorrespondence(" << bgc.pm_.left.size() << ", "
      << bgc.tm_.left.size() << ")";
  }
};



template<typename BGPlace, typename BGTransition, typename vert_t>
bool PutPlace(BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    vert_t k, BGPlace val)
{
  if (map.pm_.right.find(val)!=map.pm_.right.end()) {
    SMVLOG(BOOST_LOG_TRIVIAL(error) << "Place "<<val<<" already exists.");
    return false;
  }
  using PlaceMap=
    typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::PlaceMap;
  using VertToPlace=typename PlaceMap::left_value_type;
  typename PlaceMap::left_map::iterator inserted;
  bool success;
  std::tie(inserted, success)=map.pm_.left.insert(VertToPlace(k, val));
  if (success==false) {
    SMVLOG(BOOST_LOG_TRIVIAL(error)<<"Failed to insert a place value.");
    assert(success);
  }
  return true;
}



template<typename BGPlace, typename BGTransition, typename vert_t>
bool PutTransition(BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    vert_t k, BGTransition val)
{
  // We allow transitions which have the same key but different vertices.
  using TransitionMap=
    typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::TransitionMap;
  using VertToTrans=typename TransitionMap::left_value_type;
  typename TransitionMap::left_map::iterator inserted;
  bool success;
  std::tie(inserted, success)=map.tm_.left.insert(VertToTrans(k, val));
  if (success==false) {
    SMVLOG(BOOST_LOG_TRIVIAL(error)<<"Failed to insert a transition value.");
    assert(success);
  }
  return true;
}


/*! Get the vertex associated with a place.
 */
template<typename BGPlace, typename BGTransition, typename vert_t>
typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::vertex_descriptor
GetPvertex(const BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    BGPlace& key)
{
  auto it=map.pm_.right.find(key);
  if (it==map.pm_.right.end()) {
    SMVLOG(BOOST_LOG_TRIVIAL(error) << "Place does not exist: "<<key
      <<" map size " << map.pm_.right.size());
    return vert_t{};
  }
  return it->second;
}



/*! Get the vertex associated with a transition. There may be more
 *  than one vertex, but this returns only the first.
 */
template<typename BGPlace, typename BGTransition, typename vert_t>
typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::vertex_descriptor
GetTvertex(const BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
    BGTransition& key)
{
  auto it=map.tm_.right.find(key);
  if (it==map.tm_.right.end()) {
    SMVLOG(BOOST_LOG_TRIVIAL(error) << "Transition does not exist: "<<key);
    return vert_t{};
  }
  return it->second;
}




template<typename BGPlace, typename BGTransition, typename vert_t>
BGPlace GetPlace(const BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
  typename BiGraphCorrespondence<BGPlace,BGTransition,
      vert_t>::vertex_descriptor& key)
{
  auto it=map.pm_.left.find(key);
  if (it==map.pm_.left.end()) {
    SMVLOG(BOOST_LOG_TRIVIAL(error) << "Vertex for place does not exist: "<<key);
    return BGPlace{};
  }
  return it->second;
}




template<typename BGPlace, typename BGTransition, typename vert_t>
BGTransition GetTransition(
  const BiGraphCorrespondence<BGPlace,BGTransition,vert_t>& map,
  typename BiGraphCorrespondence<BGPlace,BGTransition,vert_t>::vertex_descriptor& key)
{
  auto it=map.tm_.left.find(key);
  if (it==map.tm_.left.end()) {
    SMVLOG(BOOST_LOG_TRIVIAL(error) << "Vertex for transition does not exist: "
      <<key);
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
  std::map<vert,std::unique_ptr<Transition>> transitions_;

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
    auto transition_vertex=add_vertex({PetriGraphColor::Transition, 0}, g_);
    bool added=PutTransition(bimap_, transition_vertex, t);

    for (auto edge : e) {
      BGPlace p=std::get<0>(edge);
      int weight=std::get<1>(edge);
      if (!this->AddEdge(transition_vertex, p, weight)) {
        added=false;
      }
    }

    transitions_.emplace(transition_vertex, std::move(transition));

    assert(added);
    return added;
  }



  bool AddEdge(vert transition, BGPlace p, int weight)
  {
    assert(g_[transition].color==PetriGraphColor::Transition);
    auto pv=GetPvertex(bimap_, p);
    assert(g_[pv].color==PetriGraphColor::Place);

    if (weight<0) {
      edge_t new_edge;
      bool success;
      std::tie(new_edge, success)=boost::add_edge(pv, transition, {weight}, g_);
      if (!success) {
        SMVLOG(BOOST_LOG_TRIVIAL(error) << "Could not add edge");
        return false;
      }
    } else {
      edge_t new_edge;
      bool success;
      std::tie(new_edge, success)=boost::add_edge(transition, pv, {weight}, g_);
      if (!success) {
        SMVLOG(BOOST_LOG_TRIVIAL(error) << "Could not add edge");
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
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<< stoch_start
        << " stochiometric coefficients start");
    auto bipartite=IsBipartitePetriGraph(g_);
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<< " is bipartite "<<bipartite);
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
      auto place_iter=bimap_.pm_.left.find(old_vertex);
      auto transition_iter=bimap_.tm_.left.find(old_vertex);
      if (place_iter!=bimap_.pm_.left.end()) {
        if (transition_iter!=bimap_.tm_.left.end()) {
          SMVLOG(BOOST_LOG_TRIVIAL(error)<<"The same vertex points both to a place "
            "and a transition.");
          assert(transition_iter==bimap_.tm_.left.end());
        }
        PutPlace(b, static_cast<int64_t>(new_vertex), place_iter->second);
      } else if (transition_iter!=bimap_.tm_.left.end()) {
        PutTransition(b, static_cast<int64_t>(new_vertex),
            transition_iter->second);

        auto trans_obj_iter=transitions_.find(transition_iter->first);
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
