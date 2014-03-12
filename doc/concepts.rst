=================
Concepts
=================


Generalized Stochastic Petri Net
---------------------------------

A generalized stochastic Petri net (GSPN) describes the relationship
among places and transitions which define a semi-Markov system.
The C++ concept of a GSPN first identifies the type of the
identifiers for places and for transitions.::

  template<typename GSPN>
  struct petri_place
  {
    typedef void type;
  };

  template<typename GSPN>
  struct petri_transition
  {
    typedef void type;
  };


Then there are two functions that express which transitions depend on
which places. The first retrieves that part of the GSPN relevant to a
particular transition, which are the input and output places, the
stochiometric coefficients associated with those edges, and the type of
token expected at each place. The second function calls a unary function
on each transition which relates to any of a given set of places.::

  template<typename GSPN>
  std::vector<std::tuple<place_t<GSPN>,int,size_t>>
  neighbors_of_transition(const GSPN& g, trans_t<GSPN> transition_id);


  auto func=[] (trans_t<GSPN> neighbor_id)->void;

  template<typename GSPN, typename Functor>
  void neighbors_of_places(const GSPN& g, const std::set<place_t<GSPN>>& place_id,
  		Functor func);


Together, these two functions determine which transitions can affect or be
affected by the firing of a given transition.

Unstated, but important, necessary conditions are that the places form
a vector space and that the transitions form a vector space. In addition,
the places connected to a transition are ordered. This is how we enforce
that they be distinguishable, unlike those in simpler Petri nets.

The last part of the GSPN concept is the transition itself, represented
only by its transition id and its distribution. The ``enabled`` function
checks the local marking to see whether a transition is enabled and,
if so, returns the distribution associated with that transition. The
Distribution object doesn't depend on the local marking.::

  template<typename GSPN, typename State, typename LocalMarking>
  std::tuple<bool,Distribution>
  enabled(const GSPN&, const State&, const LocalMarking&);

The firing of a transition is the last part of the GSPN concept.::

  template<typename GSPN, typename State, typename LocalMarking, typename RNG>
  void fire(const GSPN&, State&, LocalMarking&, RNG&);

The random number generator is necessary because firing may select
a random token from the marking at a place.

What *isn't* in this concept? There is no explicit representation of
the transitions themselves. The GSPN presents its embedded Markov matrix
as a set of distributions. Whether to construct an explicit list of
transition objects, and how to construct them, is a representational issue
for the GSPN object.



Marking
---------
The marking represents mutually-exlusive subsets of the total
state of the system so that we can represent independent 
competing processes which affect each other by modifying subsets
of the state. The marking is the set of all tokens at all places,
where the places are identified by place ids.

In this library, we permit the marking at each place to contain
only one type of token, with possible subclasses. There may be
different token types, which we call layers. The marking is
accessed in a type-safe way by specifying the token layer to retrieve.::

  template<typename... TokenContainers>
  class Marking
  {
  public:
  };

  template<typename LAYER>
  void add(Marking&, place_t<GSPN>, token_t<Marking,LAYER>);

  template<typename LAYER>
  void length(const Marking&, place_t<GSPN>);

  template<typename LAYER>
  void length(const Marking&, place_t<GSPN>, color_t<Marking,LAYER>);

  template<typename LAYER, typename UnaryOperator>
  std::tuple<std::resultof<UnaryOperator>,bool>
  get(const Marking&, place_t<Marking>, UnaryOperator& functor);

  template<typename LayerFrom, typename LayerTo>
  void
  move(const Marking&, place_t<Marking,LayerFrom>, place_t<Marking,LayerTo>, size_t);

  template<typename LayerFrom, typename LayerTo, typename UnaryOperator>
  void
  move_modify(const Marking&, place_t<Marking,LayerFrom>, place_t<Marking,LayerTo>,
  size_t, UnaryOperator& functor);


It is possible to move tokens between layers with the same token type.
A ``TokenContainer`` holds the tokens. There are two types, ``Colored``
and ``Uncolored``. The color of a token comes from a traits class::

  template<typename Token>
  token_color
  {
    typedef void type;
  }

  template<typename Token>
  unique_color
  {
    static const bool value=true;
  };

A color is unique if there will only be one token of any given color
at a place.

If we were to ask what should be stored in the marking and what in
the state, the answer is that the embedded Markov class of this library
will automatically track changes to the marking and update the 
transition distributions and enabling times, so any state whose change
necessitates a change to enabling times belongs in the marking.
Otherwise, the user will need to signal to the embedded Markov
matrix when the state has been changed by hand.

Other changes in the state of the environment, such as temperature
affects on transition rates, can be handled as distribution functions
that depend on the current semi-Markov simulation time.


GSPNState
-----------
The GSPNState has three members, the marking, which represents all tokens
at all places, the enabling time of every enabled transition, and
the current time of the semi-Markov model, which is the sum of all
transition intervals since the start of the simulation.::

  template<typename GSPN, typename Marking, typename TimeStrategy=KahanTime>
  class GSPNState
  {
  public:
    typedef Marking Marking;
    Marking marking;
    std::map<trans_t<GSPN>,double> enabling_time;
    double current_time() const;
    double add_time(double);
  };

The GSPNState re-advertises the ``Marking`` type through a public typedef.
There are two ``TimeStrategy`` classes. ``SimpleTime`` stores time in
a double, and ``KahanTime`` uses the Kahan summation algorithm to 
ameliorate roundoff error when adding small intervals to larger times.


