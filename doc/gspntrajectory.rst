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

The partial core matrix has two responsibilities. For the
underlying GSPN, it enforces the rule, central to this kind
of semi-Markov model, that distributions of transitions are
recalculated only when a transition fires which shares marking.
Sharing markings is the only permitted way to express dependence among
processes in the system. A distribution may depend arbitrarily
on external state, such as temperature, but no transition can
affect the distribution of another transition except by changing
the marking.

The partial core matrix also has a responsiblity
to present to a propagator
distributions in time of the currently-enabled competing processes.
The propagator is a function which samples these distributions to
return the next transition to fire.::

  using HazardPropagator=NonHomogeneousPoissonProcesses<TransitionKey,RandGen>;
  HazardPropagator competing_hazards;
  using GeneralPropagator=PropagateCompetingProcesses<TransitionKey,RandGen>;
  GeneralPropagator other_processes;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing_hazards, &other_processes});

It is the StochasticDynamics class that coordinates, inside, with
the partial core matrix. It will take
enabled transitions from the GSPN and give their distributions
to the appropriate propagator, depending on whether that propagator
can calculate the firing time of that distribution.
When a transition changes the marking,
the StochasticDynamics object is able to track which parts of the
marking changed and enable or disable transitions appropriately.

Finally, the main loop becomes short::

    dynamics.Initialize(&state, &rng);
    bool running=true;
    auto nothing=[](SIRState&)->void {};
    while (running) {
      running=dynamics(state);
      output_function(state);
    }
    output_function.final(state);

Here, ``output_function`` is an object which acts on the state
and gathers results.
