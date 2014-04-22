==============================================
Mathematics of the Semi-Markov Model and GSPN
==============================================



Generalized Stochastic Petri Net
-----------------------------------
There are many flavors of GSPN defined by various authors.
The following model was chosen to support complex epidemiological
simulations while still representing a semi-Markov system.
In brief, for those familiar with GSPN, this system supports
non-exponential distributions, without inhibitor arcs, with
colored tokens, without immediate transitions, permitting
simultaneous transitions.

The **what** of the GSPN is

* Places, :math:`p_1, p_2, p_3\ldots`.

* Transitions, :math:`e_1, e_2, e_3\ldots`. Each transition has
  input places, :math:`p_k\in I(e_i)` and output places,
  :math:`p_k\in J(e_i)`. For each pair of transition and input
  or output place, there is an integer stochiometric coefficient.
  There are no inhibitor places in this representation.

* Marking at each place, :math:`s_1, s_2, s_3\ldots`. The marking
  at a place is a non-negative number of tokens, :math:`q_i`. Tokens may
  have properties, :math:`\xi(q_i)`.

The **how** of the GSPN is

#. A transition becomes *enabled* when there are as many tokens on
   the input places as the stochiometric coefficients require.
   Whether the output places need to be empty is a policy choice.
   The enabling time of a transition is part of the state of the
   system.

#. An enabled transition computes from the marking on its input
   and output places, and from any properties of the tokens in that
   marking, a continuous stochastic variable which determines the
   likelihood the transition will fire at any time in the future.
   This stochastic variable need not be an exponential distribution.

#. The system samples jointly from the stochastic variables of every
   enabled transition in order to determine which transition fires
   next. This is described in grave detail later.

#. When a transition fires, it removes from the marking of its input
   places the specified number of tokens and puts them in the output
   places. Any disparity in token count destroys or creates tokens.
   Which token is moved is specified by policy to be the first,
   last, or random token in the marking of the place.

#. Transitions which do not fire fall into four classes. Those
   transitions whose input and output places are not shared with
   inputs and outputs of the fired transitions retain their
   enabling time and may not make any changes to how their
   continuous variable distributions depend on the marking.
   Any dependence between transitions must be expressed by
   sharing intput or output places.

   Those previously-enabled transitions which share input or output
   places, and for which the new marking still meets their
   enabling requirements, retain their original enabling times
   but may recalculate their stochastic variable distributions
   depending on the marking.

   Previously-enabled transitions whose input or output
   places no longer have requisite tokens for enabling become
   disabled.

   All transitions which may become enabled after the firing
   are enabled with an enabling time set to the current system
   time.



Each transition defines an independent continuous stochastic variable,
whose density is :math:`f_\alpha(\tau, t_e)` and cumulative density
is :math:`F_\alpha(\tau, t_e)`, where :math:`t_e` is the time
at which the transition was enabled. If we define a new stochastic variable,
whose value is the minimum of all of the :math:`f_\alpha`, then it
defines the time of the next transition, and the index :math:`\alpha` of
the transition that came first determines, by its transition's 
change to the marking, the new system state, :math:`j`. This means we have
a race. The racing processes are dependent only when they modify
shared marking.

If the :math:`f_\alpha` are all exponential processes, then the minimum
of the stochastic variates is analytically-tractable, and the result
is exactly Gillespie's method. Using a random number generator to sample
each process separately would, again for exponential processes, result
in the First Reaction method. The Next Reaction method is just another
way to do this sampling in a statistically-correct manner, valid only
for exponential processes.

For general distributions, the non-zero part of the core matrix,
which we will call :math:`p_{nm}(\tau)` for partial core matrix,
right after the firing of any transition, is

.. math::

   p_{nm}(\tau)=\frac{f_\alpha(\tau,t_0,t_e)}{1-F_\alpha(\tau,t_0,t_e)}\prod_\beta(1-F_\beta(s, t_0, t_e))

Here :math:`(n,m)` label states determined by how each transition, :math:`\alpha`
changes the marking. The system has, inherently, a core matrix, but we
compute its trajectory from this partial core matrix.


Semi-Markov Processes
----------------------
There are many excellent books on semi-Markov models,
notably [Howard2007]_. Our goal here is to define the generalized
stochastic Petri net as a representation of a semi-Markov
model in order to explain how the library decomposes it
into algorithms.


We define the semi-Markov process from the Markov process, following
[Pyke1961]_. A *Markov chain* is a discrete variable on a set of
states, :math:`J`. At each transition, a next state is selected
from the set :math:`J` with probability determined by the matrix
:math:`\pi_{ij}`. You are in state :math:`i` and going to one of
the states :math:`j`, with no specification of a time.

A *Markov process* introduces the time at which the next state
will be chosen. It is defined on a joint discrete and continuous
space, :math:`(J,X)`, where the continuous space runs from minus
infinity to infinity. We usually think of a Markov process as
having a continuous variable that is exponentially-distributed
in time. With a process of this sort, you can interrupt it at
any moment between transitions and still be able to say with
the same certainty, that it will likely fire within a given interval.
This property is called memorylessness. A Markov process doesn't
inherently have this restriction, though. The distribution of
its continuous stochastic variable can have any distribution,
such as Weibull, Gamma, or piecewise-continuous.

While a Markov process determines the choice of the
next state and time of the next state, a semi-Markov process
records the state of a Markov process at any time. If we track
the sum of times, :math:`S=∑Δt_{ij}`, then the semi-Markov
process tracks the state of the system at time :math:`S`.
The result is that we have a system which can have a
non-exponential sojourn in any state.

How do we simulate a trajectory from such a process?
We use the Markov process, called the *embedded* Markov process,
to select a state at each time. The *core matrix* determines the
probability of selecting the next state and time.

.. math::

   q_{ij}(\tau)=P(j, \tau|i, t_0)

There are a few choices for how to sample from a joint
discrete and continuous process. The first method is to factorize
the probability into a marginal and conditional distribution. Select 
uniformly from the marginal distribution, and then use that result
to select uniformly from the conditional distribution at that value.

.. math::

   q_{ij}(\tau)=\pi_{ij}(\infty)h_{ij}(\tau)

   q_{ij}(\tau)=\pi_{ij}(\tau)w_{i}(\tau)

The first equation marginalizes time-dependence to find the
Markov chain for the system, :math:`\pi_{ij}(\infty)`. Given a chosen
transition, the :math:`h_{ij}(\tau)` are the holding times for next
states. The second equation marginalizes over all possible next transitions
of the system to decide first when it will fire, as given by the
single waiting time, :math:`w_i(\tau)`. Then, at a particular time
it finds the probability of transtions, :math:`\pi_{ij}(\tau)`.

Imagine that we had a system with a very large state space, possibly
infinite, but that, at any point in time only a finite number of
transitions were possible. The core matrix would be infinite, but its
non-zero entries would be finite. We could simulate such a system
by looking at the non-zero entries after each transition. We would just
need to track which transitions existed at each step.




.. [Howard2007] R. A. Howard, Dynamic Probabilistic Systems: Semi-Markov and Decision Processes. Mineola, NY: Dover, 2007.

.. [Pyke1961] R. Pyke, “Markov Renewal Process: Definition and Preliminary Properties,” Ann. Math. Stat., vol. 32, no. 4, pp. 1231–1242, 1961.
