#ifndef _BREADTH_DOUBLE_BLACK_H_
#define _BREADTH_DOUBLE_BLACK_H_ 1

#include <tuple>
#include <list>
#include <map>
#include <vector>
#include "boost/log/core.hpp"
#include "logging.hpp"

namespace afidd
{
namespace smv
{

/*! Breadth first search, coloring nodes and coloring edges.
 *  In addition to nodes that are white, grey, and black, the
 *  edges can be white, grey and black, to specify when they are
 *  visited 0, 1, or 2 times. When all edges of a node are black,
 *  the node is changed to a fourth color, called Gandalf.
 *
 *           node                edge
 *  white    exists              exists
 *  grey     queued              grey parent and grey child
 *  black    children queued     grey parent and black child
 *  gandalf  all edges black     black parent and black child
 *
 *  The extra node color is useful because, during a sweep across
 *  a graph, a Gandalf node is no longer connected to any non-black
 *  node and can be removed from use for some algorithms.
 *
 *  The value_type for this Boost::iterator_range style of iterator
 *  is a tuple (int, one or two vertices in a vector). The number of
 *  vertices indicates either a node or an edge.
 */
template<typename Graph>
class BreadthDoubleBlack
{
  using Vert=typename boost::graph_traits<Graph>::vertex_descriptor;
  using VertIter=typename boost::graph_traits<Graph>::vertex_iterator;
  using Edge=typename boost::graph_traits<Graph>::edge_descriptor;
  using EdgeIter=typename boost::graph_traits<Graph>::edge_iterator;

  enum Color { White, Grey, Black, Gandalf };

public:
  using value_type=std::pair<int,std::vector<Vert>>;

private:
  const Graph& _g;
  std::list<value_type> _val; // Breadth-ordered list of verts and edges.
  std::list<Vert> _q; // queue of vertices to visit
  std::map<Vert,int> _color; // color of every vertex
  // For each vertex, number of edges which are not black (visited from
  // the vertex on each side). Starts as vertex degree.
  std::map<Vert,size_t> _degree;
  std::map<Edge,int> _edge_color; // color of every edge

public:
  BreadthDoubleBlack(const Graph& g, Vert v)
  : _g(g) {
    VertIter vb, ve;
    std::tie(vb, ve)=vertices(g);
    for (; vb!=ve; ++vb) {
      _color[*vb]=White;
      _degree[*vb]=out_degree(v, g);
    }

    EdgeIter eb, ee;
    std::tie(eb, ee)=edges(g);
    for (; eb!=ee; ++eb) {
      _edge_color[*eb]=White;
    }

    _q.push_back(v);
    _val.push_back({Grey, {v}});
    BOOST_LOG_TRIVIAL(debug) << "BFS " << Grey << " " << v;

    refill();
  }

  ~BreadthDoubleBlack() {};

  // This is a Boost::iterator_range interface. Easier than std::iterator.
  value_type& front() { return _val.front(); }

  bool empty() { return _val.empty() && _q.empty(); }

  void advance_begin(void) {
    assert(!_val.empty());
    if (_val.empty()) {
      BOOST_LOG_TRIVIAL(error) << "BFS vertices should not be empty";
      throw std::runtime_error("List of vertices should not be empty.");
    }
    _val.pop_front();
    if (_val.empty() and !_q.empty()) {
      refill();
    }
  }

  void refill() {
    Vert parent=_q.front();
    _q.pop_front();

    std::vector<Vert> verts_to_gandalf;
    using EdgeIter=typename boost::graph_traits<Graph>::out_edge_iterator;
    EdgeIter eb, ee;
    std::tie(eb, ee)=out_edges(parent, _g);
    for (; eb!=ee; ++eb) {
      auto child=target(*eb, _g);
      if (_color[child]==White) {
        _q.push_back(child);
        _color[child]=Grey;
        _val.push_back({Grey, {child}});
        BOOST_LOG_TRIVIAL(trace) << "BFS " << Grey << " " << child;
      }

      if (_edge_color[*eb]==White) {
        _edge_color[*eb]=Grey;
      }
      else if (_edge_color[*eb]==Grey) {
        _edge_color[*eb]=Black;
        _degree[parent]-=1;
        _degree[child]-=1;
        if (_degree[child]==0) {
          verts_to_gandalf.push_back(child);
        }
      }
      _val.push_back({_edge_color[*eb], {parent, child}});
      BOOST_LOG_TRIVIAL(trace) << "BFS " << _edge_color[*eb] << " "
          << parent << " " << child;
    }

    _color[parent]=Black;
    _val.push_back({Black, {parent}});
    BOOST_LOG_TRIVIAL(trace) << "BFS " << Black << " " << parent;

    for (auto& vg : verts_to_gandalf) {
      _color[vg]=Gandalf;
      _val.push_back({Gandalf, {vg}});
      BOOST_LOG_TRIVIAL(trace) << "BFS " << Gandalf << " " << vg;
    }

    if (_degree[parent]==0) {
      _val.push_back({Gandalf, {parent}});
      BOOST_LOG_TRIVIAL(trace) << "BFS " << Gandalf << " " << parent;
    }
  }
};


}
}

#endif // _BREADTH_DOUBLE_BLACK_H_
