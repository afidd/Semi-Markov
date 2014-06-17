********************
Dynamics Reference
********************

.. default-domain:: cpp

GSPNState
==============

This holds the state for the whole system.

Declaration
^^^^^^^^^^^^^^^^^^^^^^
.. class:: afidd::smv::GSPNState\<typename Mark, typename TransitionKey, typename UserState>

   The ``Mark`` is the type of the Marking class. ``TransitionKey`` is
   the transition key for the GSPN. The ``UserState`` may be omitted
   but is whatever state you would like to add to the system state.
   It will be accessible within transition as the ``UserState``
   variable and is accessible from the state as ``GSPNState.user``.

Member Functions
^^^^^^^^^^^^^^^^^
.. function:: double CurrentTime()

   Returns the current time in the system.

.. function:: void SetTime(double time)

   Sets the system time.


Members
^^^^^^^^^^^^^^^^^^

.. member:: Mark marking

   The type of this is the Mark type passed in.

.. member:: UserState user

   This is the default-instantiated value of the UserState type.

.. member:: TransitionKey last_transition

   The transition key for the last transition to fire.


Stochastic Dynamics
======================

Declaration
^^^^^^^^^^^^^^

.. class:: afidd::smv::StochasticDynamics\<typename GSPN,typename DynState,typename RandGen\>

   ``StochasticDynamics`` is responsible for advancing the state of
   the system using the given Propagators. GSPN is the type of the
   Petri net. DynState is the type of the state of the system.
   RandGen is the type of the random number generator.

Member Functions
^^^^^^^^^^^^^^^^^

.. function:: StochasticDynamics(GSPN& gspn, PropagatorVector& pv)

   The propagator vector is a type declared in this class.
   In its constructor, we pass the model for the system and the list
   of propagators which determine the next state for the system.


.. function:: Initialize(State* state, RangGen *rng)

   There may be internal arrays initialized at this stage once the model
   and state are known. This must be called before running the
   ``StochasticDynamics.``


.. function:: operator()(State& state)

   This modifies the state, in place, according to the propagators.



Continuous Propagator
=========================

The pure virtual base class for propagators of the state.
There are multiple propagators because some propagators are more
efficient than others at handling subsets of the transitions
enabled in the system.

Declaration
^^^^^^^^^^^^^

.. class:: afidd::smv::ContinuousPropagator\<typename TransitionKey, typename RNG\>

   The arguments are the transition key of the GSPN and the random number generator type.

Pure Virtual Member Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: bool Include(const TransitionDistribution\<RNG\>& distribution) const

   A propagator can examine a distribution from the core matrix to decide whether
   it will accept it.

.. function::  std::tuple\<TransitionKey, double\> Next(double now, RNG& rng) const=0

   What is the next predicted transition from this propagator.

.. function::  std::tuple\<bool,double\> Enabled(const TransitionKey& tkey) const

   Is this transition currently enabled?

.. function::  void Enable(const TransitionKey& tkey, std::unique_ptr\<TransitionDistribution\<RNG\>\>& distribution, double when, bool previously_enabled, RNG& rng)

   Enable this transition at this time.

.. function::  void Disable(const TransitionKey& tkey, double when)

   Disable this transition at this time.

.. function::  void Fire(const TransitionKey& tkey, double when, RNG& rng)

   This transition has fired.


PropagateCompetingProcesses
=============================
This subclass of ``ContinuousPropagator`` uses Gillespie's First Reaction
method ([Gillespie:1978]_) to choose the next transition. It's methods are exactly those listed
above.

Declaration
^^^^^^^^^^^^

.. class:: afidd::smv::PropagateCompetingProcesses\<typename TransitionKey, typename RNG\>

.. function:: PropagateCompetingProcesses()

   The constructor has no arguments.


NonHomogeneousPoissonProcesses
=================================
This subclass of ``ContinuousPropagator`` uses Anderson's method ([Anderson:2007]_)
to choose the next transition. It's methods are exactly those listed
above.


Declaration
^^^^^^^^^^^^

.. class:: afidd::smv::NonHomogeneousPoissonProcesses\<typename TransitionKey, typename RNG\>

.. function:: NonHomogeneousPoissonProcesses()

   The constructor has no arguments.
