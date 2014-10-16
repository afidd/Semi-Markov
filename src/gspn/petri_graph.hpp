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
#ifndef _PETRI_GRAPH_H_
#define _PETRI_GRAPH_H_ 1

#include <vector>
#include "stochnet.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/property_map/vector_property_map.hpp"
#include "marking.hpp"
#include "logging.hpp"


namespace afidd
{
namespace smv
{
  enum class PetriGraphColor : int { Unused, Place, Transition };

  struct PetriGraphVertexProperty {
    PetriGraphColor color;
    int token_layer;
  };

  struct PetriGraphEdgeProperty {
    int stochiometric_coefficient;
  };

  struct PetriGraphProperty {};

  typedef boost::adjacency_list<
    boost::listS, // VertexList container
    boost::listS, // OutEdgeList container
    boost::bidirectionalS, // directionality
    PetriGraphVertexProperty, // Vertex property
    PetriGraphEdgeProperty, // Edge property
    PetriGraphProperty, // Graph property
    boost::listS // EdgeList container
    >
    PetriBuildGraphType;




  using PetriGraphType=boost::adjacency_list<
    boost::vecS, // VertexList container
    boost::vecS, // OutEdgeList container
    boost::bidirectionalS, // directionality
    PetriGraphVertexProperty, // Vertex property
    PetriGraphEdgeProperty, // Edge property
    PetriGraphProperty, // Graph property
    boost::vecS // EdgeList container
    >;



  template<typename Graph>
  bool IsPlace(
      const typename boost::graph_traits<Graph>::vertex_descriptor& v,
      const Graph& g)
  {
    return g[v].color==PetriGraphColor::Place;
  }


  template<typename Graph>
  bool IsTransition(
      const typename boost::graph_traits<Graph>::vertex_descriptor& v,
      const Graph& g)
  {
    return g[v].color==PetriGraphColor::Transition;
  }


  template<typename Graph>
  int TokenLayer(
      const typename boost::graph_traits<Graph>::vertex_descriptor& v,
      const Graph& g)
  {
    return g[v].token_layer;
  }


  template<typename Graph>
  int StochiometricCoefficient(
      const typename boost::graph_traits<Graph>::edge_descriptor& e,
      const Graph& g)
  {
    return g[e].stochiometric_coefficient;
  }





/*! Find all places in and out of a transition, in order.

 */
template<typename Graph>
std::vector<std::tuple<int64_t,int,int>>
GraphNeighborsOfTransition(Graph& g, int64_t trans_id) {
  assert(g[trans_id].color==PetriGraphColor::Transition);
  std::vector<std::tuple<int64_t,int,int>> place_ids;

  auto initer=in_edges(trans_id, g);
  for (; initer.first!=initer.second; ++initer.first) {
    auto sc=g[*initer.first].stochiometric_coefficient;
    auto place_id=source(*initer.first, g);
    assert(g[place_id].color==PetriGraphColor::Place);
    auto level=g[place_id].token_layer;
    place_ids.push_back(std::make_tuple(place_id, level, sc));
  }
  auto oiter=out_edges(trans_id, g);
  for (; oiter.first!=oiter.second; ++oiter.first) {
    auto sc=g[*oiter.first].stochiometric_coefficient;
    auto place_id=target(*oiter.first, g);
    assert(g[place_id].color==PetriGraphColor::Place);
    auto level=g[place_id].token_layer;
    place_ids.push_back(std::make_tuple(place_id, level, sc));
  }

  return place_ids;
}

template<typename Graph>
std::vector<std::tuple<int64_t,int,int>>
GraphInputsOfTransition(Graph& g, int64_t trans_id) {
  assert(g[trans_id].color==PetriGraphColor::Transition);
  std::vector<std::tuple<int64_t,int,int>> place_ids;

  auto initer=in_edges(trans_id, g);
  for (; initer.first!=initer.second; ++initer.first) {
    auto sc=g[*initer.first].stochiometric_coefficient;
    auto place_id=source(*initer.first, g);
    assert(g[place_id].color==PetriGraphColor::Place);
    auto level=g[place_id].token_layer;
    place_ids.push_back(std::make_tuple(place_id, level, sc));
  }
  return place_ids;
}




/*! Find all transitions which depend on the Marking at a set of places.
 */
template<typename F>
void GraphNeighborsOfPlaces(PetriGraphType& g,
  const std::set<int64_t>& place_id, F func) {
  auto seen=std::set<int64_t>();

  auto report=[&seen, &func, &g] (int64_t id) {
    SMVLOG(BOOST_LOG_TRIVIAL(trace) << "transition " << id);
    assert(g[id].color==PetriGraphColor::Transition);
    if (seen.find(id)==seen.end()) {
      func(id);
      seen.insert(id);
    }
  };

  for ( auto p : place_id) {
    SMVLOG(BOOST_LOG_TRIVIAL(trace) << "neighbors of place " << p);
    assert(g[p].color==PetriGraphColor::Place);
    auto ie=in_edges(p, g);
    for (; ie.first!=ie.second; ++ie.first) {
      report(source(*ie.first, g));
    }
    auto oe=out_edges(p, g);
    for (; oe.first!=oe.second; ++oe.first) {
      report(target(*oe.first, g));
    }
  }
}




template<typename Graph>
bool IsBipartitePetriGraph(const Graph& g) {
  bool pass=true;
  using Vert=typename boost::graph_traits<Graph>::vertex_iterator;
  Vert cur_vert;
  Vert end_vert;
  std::tie(cur_vert, end_vert)=vertices(g);
  for (; cur_vert!=end_vert; ++cur_vert) {
    auto vert_color=g[*cur_vert].color;
    auto oe=out_edges(*cur_vert, g);
    for (; oe.first!=oe.second; ++oe.first) {
      if (vert_color==g[target(*oe.first, g)].color) {
        pass=false;
      }
    }
  }
  return pass;
}






/*! Count the number of stochiometric coefficients on edges.
 *  This is just one check that two graphs are not different.
 */
template<typename Graph>
size_t NumStochiometricCoefficients(const Graph& g) {
  size_t cnt{0};
  using Edge=typename boost::graph_traits<Graph>::edge_iterator;
  Edge cur_edge;
  Edge end_edge;
  std::tie(cur_edge, end_edge)=edges(g);
  for (; cur_edge!=end_edge; ++cur_edge) {
    if (g[*cur_edge].stochiometric_coefficient!=0) {
      ++cnt;
    }
  }
  return cnt;
}


} // smv
} // afidd


#endif // _PETRI_GRAPH_H_
