#ifndef _CONTACT_GRAPH_H_
#define _CONTACT_GRAPH_H_ 1


#include "boost/graph/erdos_renyi_generator.hpp"
#include "boost/graph/copy.hpp"
#include "boost/property_map/function_property_map.hpp"
#include "boost/graph/undirected_graph.hpp"


namespace afidd
{
namespace smv
{


// A contact graph needs to have an 0-based index associated with an edge
// and an id associated with a node. The id may be any time.
// The set() edge_list ensures only one edge between two nodes.
using ContactGraph=boost::adjacency_list<
                  boost::vecS, // VertexList container
                  boost::setS, // OutEdgeList container
                  boost::undirectedS, // Directed/bidirectional
                  boost::no_property, // vert prop
                  boost::no_property, // Edge property
                  boost::no_property, // Graph property
                  boost::vecS // EdgeList container
                  >;


// The edge property is there to enumerate edges.
using IndexedGraph=boost::adjacency_list<
                  boost::vecS, // VertexList container
                  boost::vecS, // OutEdgeList container
                  boost::undirectedS, // Directed/bidirectional
                  boost::no_property, // Vertex property
                  boost::property<boost::edge_index_t,size_t>, // Edge property
                  boost::no_property, // Graph property
                  boost::vecS // EdgeList container
                  >;


/*! Given a list of identifiers, make a complete contact graph.
 *  The contact graph will have a vertex property with the names.
 */
template<typename VertexName>
class CompleteGraph {
 public:
  using GraphType=boost::undirected_graph<VertexName>;

  GraphType Create(const std::vector<VertexName>& names) {
    ContactGraph graph;
    auto name_iter=names.cbegin();
    for (; name_iter!=names.cend(); ++name_iter) {
      graph.add_vertex(name_iter);
    }

    auto left=boost::vertices(n);
    for ( ; left.first!=left.second; ++left.first) {
      auto right=left.first;
      for ( ++right; right!=left.second; ++right) {
        boost::add_edge(graph, left, right);
      }
    }
    return GraphType;
  }
};


// This is used to ask boost::graph_copy not to copy the properties
// on an edge. Using just a function will not work b/c BGL looks
// for some kind of properties on it.
template<typename Graph1, typename Graph2>
struct FakeEdge
{
  template<typename Edge1, typename Edge2>
  void operator()(const Edge1& e1, Edge2& e2)
  {
  }
};


class ErdosRenyiContactGraph
{
  using ConstructType=ContactGraph;
  using ConstructVertex=boost::graph_traits<ConstructType>::vertex_descriptor;
public:
  using GraphType=IndexedGraph;
  using VertexIdProperty=boost::typed_identity_property_map<size_t>;
  using EdgeIdProperty=boost::property_map<GraphType,boost::edge_index_t>::type;

  template<typename RNGen>
  static std::tuple<GraphType,VertexIdProperty,EdgeIdProperty>
  Create(RNGen& rn_gen,
      size_t node_cnt, double edge_fraction) {
    // Use an edge list with a set in order to generate the graph without
    // having two edges connecting the same two nodes.
    using ERGen=boost::sorted_erdos_renyi_iterator<RNGen,ConstructType>;

    BOOST_LOG_TRIVIAL(debug) <<
        "ErdosRenyiContactGraph::create ergen edge_fraction="
        << edge_fraction;
    bool self_loops=false;
    ConstructType graph(ERGen(rn_gen, node_cnt, edge_fraction, self_loops),
        ERGen(), node_cnt);

    BOOST_LOG_TRIVIAL(debug)
        << "ErdosRenyiContactGraph::create denumerate";
    // Then copy it to a graph where the vertex descriptor is an ordinal.
    std::map<ConstructVertex const,size_t> denumerate;
    boost::associative_property_map<std::map<ConstructVertex const,size_t>>
        denumerate_map(denumerate);
    size_t cnt=0;
    
    for (auto vi=boost::vertices(graph); vi.first!=vi.second; vi.first++) {
      denumerate[*vi.first]=cnt++;
    }

    FakeEdge<ConstructType,GraphType> fake_e;

    GraphType indexed_copy;
    boost::copy_graph<ConstructType,GraphType>(graph, indexed_copy,
      boost::vertex_index_map(denumerate_map).edge_copy(fake_e));

    VertexIdProperty vid;

    size_t edge_idx=0;
    EdgeIdProperty edge_idx_map=boost::get(boost::edge_index_t(), indexed_copy);
    for (auto ei=boost::edges(indexed_copy); ei.first!=ei.second; ++ei.first) {
      boost::put(edge_idx_map, *ei.first, edge_idx);
      ++edge_idx;
    }

    return std::make_tuple(indexed_copy, vid, edge_idx_map);
  }
};

}
}



#endif
