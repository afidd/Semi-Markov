========================
ExplicitTransitions
========================


The ExplicitTransitions class implements the GSPN concept in two parts.
The Petri net, with stochiometry, is represented by a graph from
the Boost Graph Library. The transitions of the GSPN are represented
by Transition objects, each of which can calculate whether a transition
is enabled and what it does when it fires.

There are three steps to making a GSPN using ExplicitTransitions:
Build the graph, make a marking, then add Transitions and the graph
to the ExplicitTransitions object, which represents the GSPN concept.

The Boost graph models the bipartite Petri net by having vertices
with two colors, Place and Transition. It is a bidirectional graph,
meaning it is directed but with O(1) access to the inputs of a transition,
which a simple directional graph would not offer. The input and output
edges of all transitions are ordered, because a transition in a GSPN
can depend arbitrarily on the marking. This means a transition can
as for its second input instead of asking for the input of a particular
place. The Boost graph implementation is therefore::

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

    struct PetriGraphProperty {};

    using PetriGraphType=boost::adjacency_list<
        boost::vecS, // VertexList container
        boost::vecS, // OutEdgeList container
        boost::bidirectionalS, // directionality
        PetriGraphVertexProperty, // Vertex property
        PetriGraphEdgeProperty, // Edge property
        PetriGraphProperty, // Graph property
        boost::vecS // EdgeList container
        >;

A single graph can have different types of tokens. Maybe one represents
available food and another represents a consumer of that food. Each
token type is stored in a separate layer within the marking. The
token_layer vertex property is a sanity check that the correct token
is at the correct place.
The actual GSPN implementation stores a copy of this graph, not a reference.

Decide what types this model uses. We can do this first to get the dependencies among types out of the way.::

   // Everybody always needs a random number generator.
   using RandGen=boost::random::mt19937;
   // The Net.
   using ModelGraph=PetriGraphType;
   // The type of token.
   struct IndividualToken {};
   // The marking on the net.
   using ModelMarking=Marking<place_t<ModelGraph>, Uncolored<IndividualToken>>;
   // The state of the system.
   using ModelState=GSPNState<ModelGraph, ModelMarking>;
   // Now we can make the GSPN type.
   using ModelGSPN=ExplicitTransitions<LocalMarking<Mark>,
     ModelState,ModelGraph,RandGen>;

Each transition is derived from a base class, ExplicitTransitions::Transition.
Different transitions in the model are different derived classes of this.
A transition defines two methods. The enabled() method returns whether a
transition is enabled and what its distribution is if it is enabled.
The fire() method determines how to move and change tokens when they fire.
Both read and/or modify the current system state and the local marking, which
has convenience functions to read or move tokens automatically given
stochiometry of the transition.::

   class ModelMovement : public ModelGSPN::Transition
   {
      virtual std::pair<bool, std::unique_ptr<Dist>>
      enabled(const SIRState& s, const LocalMarking<Mark>& lm) const override
      {
        if (lm.template input_tokens_sufficient<0>())
        {
          return {true, std::unique_ptr<ExpDist>(new ExpDist(1.0))};
        }
        else
        {
          return {false, std::unique_ptr<NoDist>(new NoDist())};
        }
      }

      virtual void fire(SIRState& s, LocalMarking<Mark>& lm,
          RandGen& rng) const override
      {
        BOOST_LOG_TRIVIAL(debug) << "Fire movement " << lm;
        lm.template transfer_by_stochiometric_coefficient<0>(rng);
      }
   };

There are likely many types of transitions in any given GSPN.

Finally the types are decided. Now we create the objects.
The order of instantiation depends on the problem For this exercise,
assume there are lists containing places, transitions, and edges.::

    // Initialize the Boost Graph with a number of vertices.
    ModelGraph graph(37);
    ModelGSPN gspn(graph); // Make the GSPN.

    // Identify each place in the graph.
    PetriGraphVertexProperty vprop;
    for (auto place_idx : places)
    {
      vprop.color=PetriGraphColor::Place;
      vprop.token_layer=0;
      gspn.graph[place_idx]=vprop;
    }

    // Identify each transition.
    for (auto transition_idx : transitions)
    {
      vprop.color=PetriGraph::Transition;
      gspn.graph[transition_idx]=vprop;
      gspn.transitions.emplace(transition_idx,
        std::move(std::unique_ptr<ModelGSPN::Transition>(
        new ModelMovement())));
    }

    // Edges get stochiometry.
    for (auto edge : edges)
    {
      auto left_vertex=std::get<0>(edge);
      auto right_vertex=std::get<1>(edge);
      auto stochiometry=std::get<2>(edge);
      add_edge(left_vertex, right_vertex, stochiometry, gspn.graph);
    }

By the end, a GSPN with places and transitions is complete.
