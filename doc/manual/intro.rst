==================================================
Introduction
==================================================

The ``Semi-Markov`` library to simulates trajectories of semi-Markov
processes in continuous time. Its main use is to specify and simulate
represents these processes as Generalized Stochastic Petri Nets (GSPN).
These nets are used in many areas, from reliability analysis and
computer performance modeling to chemical systems. While it would
support those uses, as well, the flavor of GSPN chosen for this
library is tuned to model epidemiological and ecological models.
It can model the stages of a disease, life history of an individual,
movement of individuals, or all of these at once.

There are a number of common techniques for epidemiological
and ecological dynamical modeling, such as

* deterministic differential equations

* stochastic differential equations, and

* individual-level models in discrete time.

The GSPN can be a specification for either a discrete time system or,
as used here, continuous time system. We emphasize that, given
field measurements at the individual level, a GSPN can produce
a model most faithful to the measured data, so that it serves as
specification even when solution of larger systems is not tractable.
At the same time, this library uses algorithms to make possible
enormous stochastic, continuous-time simulations.

This library was created the Analytical Framework for Infectious Disease
Dynamics (AFIDD) group at Cornell University in conjunction with
the USDA. The work was funded by the Department of Homeland
Security, under the terms of which contract it is placed into the Public
Domain.

What is a GSPN?
-----------------

A Generalized Stochastic Petri Net (GSPN) is a formal way to specify a
a system of interacting, competing processes. Different organisms
can compete, but, for this system, the likelihood of infecting
a neighbor versus the likelihood of recovery are seen as competing,
as well.

Define a system by placing *tokens* at *places,* the way you would
put checkers on a game board. Each place represents a sub-state of
the system, such as herd of animals. Five tokens on a place representing
a herd means the herd has five animals.

*Transitions* compete to move the tokens. Each transition is
an independent process. (We explain later how and why independent processes
are able to represent biological processes that are clearly dependent.)
Only transitions change the state. Each one triggers according to its
own internal clock. This library can model non-exponential distributions
of firing times.

To get a flavor, try the following small examples.


Frogs on Lily Pads
---------------------

Susceptible-Infected-Susceptible
----------------------------------


Dairy Herd Example
-------------------

As an example, let's look at a simple model for the farm management
of dairy cows, following [viet2004]_.
There are four groups, calves, heifers before breeding
(heifer 1), heifers ready for breeding (heifer 2), and dairy cows
which give birth to new calves. We might make a sketch as shown here.

.. image:: images/bvd_gist.svg
   :scale: 50%
   :alt: Cow management has four stages, calf, heifer 1, heifer 2, and dairy.
   :align: center

We may have several different goals for this model. We may want to ask
how quickly a disease might spread through a herd, on average. We may
want to parameterize a differential equation model for changes in
herd sizes given economic conditions. The data for this model, though,
come in the form of charts of how many days it took a particular set
of heifers to give birth after their first insemination. It comes in
the form of rules that farmers with too many calves to fit in the pen
sell the rest. We are going to use a GSPN to express this model
in terms of the measured quantities.

The chart above shows a set of states for the cow, but it isn't clear,
for instance, about the number of ways a cow can leave the farm.
There can be two ways a cow can *transition* to leaving the
farm: sale or death.





.. [viet2004] A.-F. Viet, C. Fourichon, H. Seegers, C. Jacob, and C. Guihenneuc-Jouyaux, “A model of the spread of the bovine viral-diarrhoea virus within a dairy herd.,” Prev. Vet. Med., vol. 63, no. 3–4, pp. 211–36, May 2004.
