==============================================
Mathematics of the Semi-Markov Model and GSPN
==============================================

There are many excellent books on semi-Markov models,
notably [Howard2007]_. Our goal here is to define the generalized
stochastic Petri net as a representation of a semi-Markov
model in order to explain how the library decomposes it
into algorithms.

Semi-Markov Processes
----------------------

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


GSPN
--------
The generalized stochastic Petri net, or GSPN, specifies a semi-Markov
system not by specifying the core matrix but by specifying long-lived,
competing processes from which, after any transition,
we can determine the non-zero part of the semi-Markov core matrix.

For the semi-Markov system, the state was chosen from a set :math:`J` and
enumerated with :math:`i` and :math:`j`. For the GSPN, the state is
partitioned into mutually-exclusive sub-states, each identified with
a *place.* The state of the system is the number of *tokens* at each
places, along with any state within the tokens. The whole state
is called the *marking.* A single state :math:`i` is now associated with
the number and values of tokens at every places in the marking.

A set of *transitions* defines the next possible states of the system.
Each transition can move and modify tokens in the marking, so the set
of tokens moved associates a particular transition with a particular
change of state. Each transition represents a competing process
which is independent over the life of the transition, meaning
it is a continuous stochastic variable whose distribution is
independent of any change in the marking which does not disable
the transition. (It may, however, change due to other changes of
state, such as inhomogenous factors. Think of dependence on weather.)
A transition is enabled by the presence and/or value
of tokens at places; it depends on a subset of the total marking.
It can be disabled when other transitions move or modify tokens. This
is how we model dependence among physical processes when using
transitions whose distributions are independent between the firing
of transitions.

Because each transition depends upon a subset of the places, we represent
the relationship as a bipartite graph, :math:`(P,T)`. Many models,
such as chemical master equations, move tokens around in fixed numbers
for each transition (or reaction), according to stochiometry,
so we label each edge with the number of tokens a transition consumes
and a number it produces, as stochiometric coefficients. A transition
is then enabled when there are sufficient input tokens.

We were able to sample from a semi-Markov model by factorizing
the core matrix. How do we represent, and sample from, the GSPN
representation?

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

   p_{nm}(\tau)=\frac{f_\alpha(\tau,t_0,t_e)}{1-F_\alpha(\tau,t_0,t_e)}\prod_\beta\int_{t_0}^∞ 1-F_\beta(s, t_0, t_e)ds

Here :math:`(n,m)` label states determined by how each transition, :math:`\alpha`
changes the marking. The system has, inherently, a core matrix, but we
compute its trajectory from this partial core matrix.

We consider the GSPN model to be:

* The set of places.
* The set of transitions.
* For each transition,
   * The set of places on whose marking the enabling
     and distribution depends.
   * The distribution itself, to be calculated at any time :math:`t_0`.
   * How the transition modifies tokens when it fires, and the places
   	 whose markings change as a result.

The state of the system is:

* The tokens at each place in the marking.
* The enabling time of every transition.

These provide a roadmap for library implementation.



.. [Howard2007] R. A. Howard, Dynamic Probabilistic Systems: Semi-Markov and Decision Processes. Mineola, NY: Dover, 2007.

.. [Pyke1961] R. Pyke, “Markov Renewal Process: Definition and Preliminary Properties,” Ann. Math. Stat., vol. 32, no. 4, pp. 1231–1242, 1961.
