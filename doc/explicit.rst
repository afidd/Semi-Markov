===============================================
Using ExplicitTransitions to Construct a GSPN
===============================================


The ExplicitTransitions class implements the GSPN concept.

The Type of the GSPN Model
------------------------------

The author of the code supplies building blocks to this class:

* A token type. It should have an ostream ``operator<<`` defined.
  It may have any necessary data, such as a birth date
  or size of truck. It must have a default constructor. The simplest
  choice is to construct a `POD type <http://en.wikipedia.org/wiki/C++11#Modification_to_the_definition_of_plain_old_data>`_.

* A storage container for tokens at a place. The allowed option
  is Uncolored<TokenType>, which is an ``std::vector.``

* A PlaceKey, which is another POD type to uniquely identify
  a place in the GSPN. It needs less-than, equal-to, and print
  operators.

* A TransitionKey, which uniquely identifies transitions in the
  GSPN. It needs less-than, equal-to, and print
  operators.

* UserState. The system state consists of marking, enabled
  transitions, and the absolute time. The optional UserState
  struct will add arbitrary members to the state for convenience.
  If a transition were to modify this state such that other
  transitions would have new different enabling or distributions,
  the system would no longer be semi-Markov. This state is
  a good place to put constant parameters or links to
  external data, such as weather.

* A random number generator. The system assumes a single
  generator for everything. This class supports the C++
  ``random`` requirements.

An example of a ``PlaceKey`` from an SIR model shows that defining
such a struct is verbose but painless::
   
    struct SIRPlace
    {
      size_t disease;
      size_t individual;

      SIRPlace()=default;
      SIRPlace(size_t d, size_t i)
      : disease(d), individual(i)
      {}

      friend inline
      bool operator<(const SIRPlace& a, const SIRPlace& b)
      {
        return lazy_less(a.disease, b.disease, a.individual,
          b.individual);
      }

      friend inline
      bool operator==(const SIRPlace& a, const SIRPlace& b)
      {
        return (a.disease==b.disease)&& (a.individual==b.individual);
      }

      friend inline
      std::ostream&
      operator<<(std::ostream& os, const SIRPlace& cp)
      {
        return os << '(' << cp.disease << ", " << cp.individual<<')';
      }
    };


The template parameters for ``ExplicitTransitions`` are as above::

  typedef LocalMarking<Uncolored<TokenType>> Local;
  typedef ExplicitTransitions<PlaceKey, TransitionKey, Local,
                              RandGen, UserState> GSPN;

The ``LocalMarking`` class is the marking of input and output tokens
for a particular transition. Instead of using ``PlaceKey`` to identify
particular sets of tokens, it assumes the transition can access its
edges with an ordinal (first, second, third). The idea is to define
transitions which don't need to know how their marking is indexed and
stored.

Define Transitions
-------------------

A transition has three responsibilities. Determine when it
is enabled, choose a continuous variable distribution function
for its firing time, and modify the marking when it fires.
Enabling and calculating a distribution happen at the same time,
so they are in one function.

A transition inherits from the ``ExplicitTransition`` template.
It acts on the local marking to produce distributions or move
tokens::

    using Transition=ExplicitTransition<Local,RandGen,WithParams>;
    using Dist=TransitionDistribution<RandGen>;

    class InfectNeighbor : public Transition
    {
      virtual std::pair<bool, std::unique_ptr<Dist>>
      enabled(const UserState& s, const Local& local_marking,
              double te, double t0) const override
      {
        if (local_marking.template input_tokens_sufficient<0>())
        {
          return {true, std::unique_ptr<Dist>(
              new WeibullDistribution(s.params.at(0), te))};
        }
        else
        {
          return {false, std::unique_ptr<Dist>(nullptr)};
        }
      }

      virtual void fire(UserState& s, Local& local_marking,
                        RandGen& rng) const override
      {
        local_marking.template transfer_by_stochiometric_coefficient<0>(rng);
      }
    };

The ``enabled()`` method's parameters are

* **UserState** - This is the same as specified above. It could
  include parameters or a pointer to inhomogeneous drivers of the system.

* **local_marking** - Most of the work here is manipulation of the
  local marking. It has methods to add a token, remove a token, move
  a token, read a value from a token, or modify a token. It also has
  convenience methods to move all tokens associated with a marking.
  The local_marking contains stochiometric coefficients for the
  transition.

* **enabling_time** - If the transition was previously-enabled,
  this is the previous enabling time. Otherwise, it is the current
  absolute time of the GSPN.

* **current time** - The current absolute time of the system.

For a newly-enabled transition, the current time and enabling time
will be the same.



Using ExplicitTransitionsBuilder
----------------------------------
The implementation of ``ExplicitTransitions`` requires that
we construct it with a builder which then produces the GSPN
object. The builder has just three methods, `add_place()`,
`add_transition()`, and `build()`. It checks that transitions
have, as inputs and outputs, places which exist. It ensures
every PlaceKey and TransitionKey is unique.

For example::

    GSPN
    build_system(size_t individual_cnt)
    {
      BuildGraph<GSPN> bg;
      using Edge=BuildGraph<GSPN>::PlaceEdge;

      enum { s, i, r };

      for (size_t ind_idx=0; ind_idx<individual_cnt; ind_idx++)
      {
        for (size_t place : std::vector<int>{s, i, r})
        {
          bg.add_place({place, ind_idx}, 0);
        }
      }

      for (size_t left_idx=0; left_idx<individual_cnt-1; left_idx++)
      {
        bg.add_transition({left_idx, left_idx, 0},
          {Edge{{i, left_idx}, -1}, Edge{{r, left_idx}, 1}},
          std::unique_ptr<SIRTransition>(new Recover())
          );

        for (size_t right_idx=left_idx+1; right_idx<individual_cnt; right_idx++)
        {
          SIRPlace left{i, left_idx};
          SIRPlace rights{s, right_idx};
          SIRPlace righti{i, right_idx};

          bg.add_transition({left_idx, right_idx, 0},
            {Edge{left, -1}, Edge{rights, -1}, Edge{left, 1}, Edge{righti, 1}},
            std::unique_ptr<SIRTransition>(new InfectNeighbor()));

          SIRPlace lefts{s, left_idx};
          SIRPlace lefti{i, left_idx};
          SIRPlace right{i, right_idx};

          bg.add_transition({right_idx, left_idx, 0},
            {Edge{right, -1}, Edge{lefts, -1}, Edge{right, 1}, Edge{lefti, 1}},
            std::unique_ptr<SIRTransition>(new InfectNeighbor()));
        }
      }

      // std::move the transitions because they contain unique_ptr.
      return std::move(bg.build());
    }

Create Marking and State
----------------------------------
The last step is to create the marking and state.
There is one implementation of a marking in the library.
While we defined a `PlaceKey` above, the `ExplicitTransitions`
class uses a different place and transition key internally,
chosen by the Boost Graph Library implementation. Both are
just type `size_t`. The marking and state are therefore::

  using Mark=Marking<size_t, Uncolored<IndividualToken>>;
  using State=GSPNState<Mark,UserState>;

  State state;

If the `UserState` had a member named `parameters`, then we could access
it as `state.user.parameters`.

How do we initialize the marking? Unfortunately, the marking doesn't
use our `PlaceKey`, so we have to add a translation step from
the `PlaceKey` we know to the `size_t` we don't::

  size_t susceptible=gspn.place_vertex(
      PlaceKey{Disease::Susceptible, individual_idx});
  add<0>(state.marking, susceptible, IndividualToken{});

The GSPN remembers the original `PlaceKey` and will translate
for us. The second line adds a new token to the marking.

That's everything that defines the model and the state of the system.
We made places, transitions, and a marking. The next step is
to simulate a trajectory of the system.

