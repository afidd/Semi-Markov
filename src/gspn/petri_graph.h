#ifndef _PETRI_GRAPH_H_
#define _PETRI_GRAPH_H_ 1

#include <vector>
#include "stochnet.h"
#include "boost/log/core.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/property_map/vector_property_map.hpp"
#include "marking.h"
#include "logging.h"


namespace afidd
{
namespace smv
{
  enum class PetriGraphColor : int { Unused, Place, Transition };
  
  struct PetriGraphVertexProperty
  {
    PetriGraphColor color;
    size_t token_layer;
  };

  struct PetriGraphEdgeProperty
  {
    int stochiometric_coefficient;
  };

  struct PetriGraphProperty
  {
  };

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
  bool is_place(
      const typename boost::graph_traits<Graph>::vertex_descriptor& v,
      const Graph& g)
  {
    return g[v].color==PetriGraphColor::Place;
  }


  template<typename Graph>
  bool is_transition(
      const typename boost::graph_traits<Graph>::vertex_descriptor& v,
      const Graph& g)
  {
    return g[v].color==PetriGraphColor::Transition;
  }


  template<typename Graph>
  size_t token_layer(
      const typename boost::graph_traits<Graph>::vertex_descriptor& v,
      const Graph& g)
  {
    return g[v].token_layer;
  }


  template<typename Graph>
  size_t stochiometric_coefficient(
      const typename boost::graph_traits<Graph>::edge_descriptor& e,
      const Graph& g)
  {
    return g[e].stochiometric_coefficient;
  }





/*! Find all places in and out of a transition, in order.

 */
std::vector<std::tuple<size_t,size_t,int>>
neighbors_of_transition(PetriGraphType& g, size_t trans_id)
{
  assert(g[trans_id].color==PetriGraphColor::Transition);
  std::vector<std::tuple<size_t,size_t,int>> place_ids;

  auto initer=in_edges(trans_id, g);
  for (; initer.first!=initer.second; ++initer.first)
  {
    auto sc=g[*initer.first].stochiometric_coefficient;
    auto place_id=source(*initer.first, g);
    assert(g[place_id].color==PetriGraphColor::Place);
    auto level=g[place_id].token_layer;
    place_ids.push_back(std::make_tuple(place_id, level, sc));
  }
  auto oiter=out_edges(trans_id, g);
  for (; oiter.first!=oiter.second; ++oiter.first)
  {
    auto sc=g[*oiter.first].stochiometric_coefficient;
    auto place_id=target(*oiter.first, g);
    assert(g[place_id].color==PetriGraphColor::Place);
    auto level=g[place_id].token_layer;
    place_ids.push_back(std::make_tuple(place_id, level, sc));
  }

  return place_ids;
}





/*! Find all transitions which depend on the Marking at a set of places.
 */
template<typename F>
void neighbors_of_places(PetriGraphType& g,
  const std::set<size_t>& place_id, F func)
{
  auto seen=std::set<size_t>();

  auto report=[&seen, &func, &g] (size_t id)
  {
    BOOST_LOG_TRIVIAL(trace) << "transition " << id;
    assert(g[id].color==PetriGraphColor::Transition);
    if (seen.find(id)==seen.end())
    {
      func(id);
      seen.insert(id);
    }
  };

  for ( auto p : place_id)
  {
    BOOST_LOG_TRIVIAL(trace) << "neighbors of place " << p;
    assert(g[p].color==PetriGraphColor::Place);
    auto ie=in_edges(p, g);
    for (; ie.first!=ie.second; ++ie.first)
    {
      report(source(*ie.first, g));
    }
    auto oe=out_edges(p, g);
    for (; oe.first!=oe.second; ++oe.first)
    {
      report(target(*oe.first, g));
    }
  }
}

} // smv
} // afidd


#endif // _PETRI_GRAPH_H_
