#ifndef _GSPN_H_
#define _GSPN_H_ 1



namespace afidd
{
namespace smv
{

namespace detail
{
  struct NoExtraState {};
}


template<typename PetriNet>
struct petri_place
{
  typedef void type;
};


template<typename PetriNet>
struct petri_transition
{
  typedef void type;
};


template<typename PetriNet>
struct petri_graph
{
  typedef void type;
};



// These are for shorthand.
template<typename PetriNet>
using trans_t=typename petri_transition<PetriNet>::type;

template<typename PetriNet>
using place_t=typename petri_place<PetriNet>::type;

template<typename PetriNet>
using graph_t=typename petri_graph<PetriNet>::type;



// These are part of the concept.

/*! Find all places in and out of a transition, in order.
template<typename PN>
std::vector<place_t<PN>>
neighbors_of_transition(graph_t<PN>& g, trans_t<PN> trans_id)
*/

/*! Find all transitions which depend on the Marking at a set of places.
template<typename PN, typename F>
void neighbors_of_places(graph_t<PN>& g,
  const std::set<place_t<PN>>& place_id, F func)
*/

/* And it returns a transition when given a transition id. */

}
}




#endif //_GSPN_H_
