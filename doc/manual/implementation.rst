*********************
Implementation
*********************

Design
=======

The two most important points about implementation of this
library are that the math matters and that we weren't kidding
about the layers.

1. The ``ExplicitTransitions`` class represents a GSPN as an
   explicit bipartite graph of places and transitions. It would
   be possible to create other classes to represent a GSPN, for instance,
   one which constructs it on-the-fly as the marking changes.

2. That representation conforms to a GSPN API, which requires
   that, given a place, the GSPN can say which transitions are affected
   by a change to that place and, given a transition, the GSPN
   can say which are input and output places.

3. The ``PartialCoreMatrix`` constructs the nonzero part of
   :math:`q_{ij}(\tau)` at every
   time step, which means it responds to the firing of a transition
   by asking which transitions have become enabled, disabled, or
   modified by the change of state due to the last transition's firing.
   It is the requirements of the ``PartialCoreMatrix`` that determine
   the GSPN API.

4. The list of enabled transitions, presented by the ``PartialCoreMatrix,``
   is the Semi-Markov API.

5. The ``StochasticDynamics`` class examines the list of enabled transitions
   in order to select the next one. It uses Propagator classes to
   perform the Gillespie First Reaction method and
   Anderson's method [Anderson:2007]_. We may implement more of these
   because different dynamics can perform better on different problems.
   Each dynamics has two responsibilities, to cache the list of enabled
   distributions and to choose among them for the next to fire. It collaborates
   with the ``PartialCoreMatrix`` which adds, removes, or updates transition
   distributions from the cache.

Separation of responsibility, induced by identification of the mathematical
models, reduces complexity in the layers. The various dynamics, for instance,
are remarkably simple, as is the explicit representation of the
GSPN. Any bodies are buried in the ``PartialCoreMatrix`` because it contains
logic for two sensitive problems, how the list of transitions is kept
current between firings and how we decide which transitions could possibly
be affected by a particular change of state. For instance, if a transition
was enabled, and it will remain enabled, it may still choose to
have a new distribution of firing times due to changes of state.


Coding Style
===============

The C++ is C++11, with lots of templates
and lots of lambda functions. Scott Meyers has a nice tutorial on
the subject [Meyers:2013]_. The language is unfortunately complicated.

The ``main()`` for a program that uses the Semi-Markov library will
have three parts. First it will construct types with either ``typedef``
or ``using declarations`` which create concrete types from the
templated ones. Then it will create instances of those types,
and finally there is a short main loop.

The types have a dependency graph, which is drawn here with arrows
going in the opposite direction:

.. figure:: images/type_dependency.*
   :width: 300px

This means you instantiate a ``PlaceKey`` and ``TransitionKey`` before
an ``ExplicitTransition.`` The examples show some choice of the order
for type instantiation. Once this is done, there are few dependencies
for instantiation of objects. For instance, the initializer for
the ``StochasticDynamics`` needs to know the GSPN and the propagators.

