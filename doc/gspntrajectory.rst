=================================================
Simulating the Trajectory of a Semi-Markov GSPN
=================================================

Given a representation of a GSPN, what steps will simulate
a trajectory from that GSPN? As a semi-Markov model, the GSPN
tracks in time the state of a system and defines the possible
transitions to next states of the system. The core matrix of a
semi-Markov model is just that, the probability of any next state
and the probability of when that a particular state happens.
We simulate a trajectory by looking at only the transitions which
are enabled for the current state. This is the *partial core matrix.*
Given this partial core matrix, a propagator selects the next state.

The PartialCoreMatrix class type is constructed from the
GSPN type, the definition of the state, and our chosen
random number generator.::

  template<typename GSPN, typename State, typename Rand>
  class PartialCoreMatrix;

Given a GSPN, its responsiblity is to present to a propagator
distributions in time of the currently-enabled competing processes.
The propagator is a function which samples these distributions to
return the next transition to fire.::

  auto next=propagate_competing_processes(partial_core, input_string, rng);

We run this repeatedly in our main loop. The input_string is a
functor passed into the propagator which can change the state of
the system, or do nothing at all. It can be used to set initial
system state, change parameters of the system, or anything else.
When this functor changes the state of the marking, in particular,
the PartialCoreMatrix object is able to track which parts of the
marking changed and enable or disable transitions appropriately.

As a whole, sampling the trajectory of a GSPN might look like this::

  PartialCoreMatrix<GSPN,State,Rand> core(gspn, state);
  auto initial_case=[&from_place, &to_place](State& state)->void {
    // Moves a token from one place to another.
    move<0,0>(state.marking, from_place, to_place, 1);
  }
  auto next=propagate_competing_processes(core, input_string, rng);

  auto nothing=[](State& state)->void {};

  for ( ;
      std::get<1>(next)<std::numeric_limits<double>::max();
      next=propagate_competing_processes(core, nothing, rng)
      )
  {
  	std::cout << "transition " << std::get<0>(next) << " at time "
  	  << std::get<1>(next) << std::endl;
  }

We have not described here observations on the state.
