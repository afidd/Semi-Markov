
*****************
Background
*****************

This addresses two main points, how to specify a model for
the library using distributions defined by hazards and
why such a specification, with its initial conditions,
is sufficient to define the trajectory for a model.


The Hazard from Survival Analysis
==================================

Discrete case
----------------

The discrete case is much easier to understand than the continuous
case because it can be explained without employing any results from
calculus.  Throughout this section, :math:`\bf{X}` will be assumed to
real-valued random variable.  For example, :math:`\bf{X}` could
represent latent periods for measles.

It frequently happens that random samples of the real valued variables
such as :math:`\bf{X}` are actually analyzed on a discrete scale.
For example Stocks' data on latent periods of measles in
:ref:`latent_period` is based on daily visits by patients.  

The (cumulative) distribution of :math:`\bf{X}` is defined as

.. math:: F_{X}(k) = \mathcal{P}[x \le k]

assuming :math:`F_{X}(\infty) = 1`.  The 
density can be expressed as the difference in adjacent values of the 
distribution

.. math:: 
   :nowrap:

   \begin{eqnarray}
   f_{X}(k) & = & \mathcal{P}[X=k] \\
            & = & \mathcal{P}[X\le k] - \mathcal{P}[X \le k-1 ] \\
	    & = & F_{X}(k) - F_{X}(k-1)
   \end{eqnarray}

For Stocks' data in :ref:`latent_period`, the density at day :math:`k`
should be interpreted as the probability of the appearance of symptoms
since the previous visit on day :math:`k-1`.

The *hazard* is defined as the conditional probability that the value
of a random selection from :math:`\bf{X}` equals :math:`k` given it
this value is already known to exceed :math:`k-1`.  Using the usual
rules for computing conditional probabilities, the hazard is given by
the following ratio

.. math:: 

   \begin{eqnarray}
   h_{X}(k) & = & \mathcal{P}[X=k\; |\; k-1<X] \\
            & = & {\frac{f_{X}(k)}{1 - F_{X}(k-1)}}
   \end{eqnarray}

In the case of Stocks' data, the hazards shown in
:ref:`latent_period_hazard` would correspond to the probability of
symptoms appearing at day :math:`k` given that the patient had not
displayed symptoms at any previous visit.  As time goes on, patients
who have already developed symptoms effectively reduce the pool of
patients in the study who are still in a state where they might first
present symptoms on day :math:`k`.  This is the origin of the term in
the denominator.

.. _latent_period_hazard:

.. figure:: images/Stocks_hazard.*
   :width: 500 px
   :align: center

   Figure 2.  Estimated hazards of latent periods for measles in
   London circa 1931

On any given day, the hazard for latent periods can be interpreted as
the rate of appearance of symptoms per asymptomatic (infected but not
yet symptomatic) patient per day.  For example, the hazard inferred
from the Weibull distribution is approximately :math:`0.15` on day 10.
In other words, 15% of the patients that are asymptomatic on day 9
will present symptoms when examined on day 10.

.. figure:: images/stocks_person.*
   :width: 200 px
   :align: center

   Figure 3. Each participant of the Stocks study could either
   become symptomatic or leave the study. Focusing on the
   hazard accounts for the effect of those who leave.

This interpretation is extremely important because it connects a
hazard with a rate for a specific process, and that rate has well
defined units of measurement.  In addition, it clarifies how rate
parameters should be estimated from observational data.  Failure to
account for the shrinking pool over time is commonplace.  In this case
it would lead to a systematic errors in the estimation of process
rates, especially at long times when the depletion effect is most
pronounced.



Continuous case
--------------------

The random variable :math:`\bf{X}` is again assumed to be a
real-valued, but the measurements will not be binned as above.
The cumulative distribution not an integer :math:`k` but a continuous
time interval, :math:`\tau`.

.. math::

   F_X(\tau)=P[x\le\tau]

assuming :math:`F_X(\infty)=1`. The density is the derivative
of the cumulative distribution. The concept of the hazard is
part of survival analysis, where survival is
:math:`G_X(\tau)=1-F_X(\tau)`, and represents the probability
the random variable, a time interval, is longer than :math:`\tau`.
One expression for the hazard is that the density of the random
variable is equal to the probability it survives to a time :math:`\tau`
multiplied by the hazard rate for firing at time :math:`\tau`, or, in
probabilities,

.. math::

   P[\tau<x\le\tau+d\tau]d\tau=P[\tau<x]P[\tau<x\le\tau+d\tau+d\tau|\tau<x].

Writing this same equation with its almost-sure equivalents defines
the continuous hazard, :math:`\lambda_X(\tau)`,

.. math::
  
   f_X(\tau)=G_X(\tau)\lambda_X(\tau).

This is a rearrangement away from the definition of the discrete case.


Finite State Machines Generate Trajectories
============================================

This library accepts a specification of a model in terms of
hazards, an initial condition, and produces trajectories.
This set of high-level steps to simulation (specify, initialize,
step) has a well-defined abstraction called a *finite state machine.*
It isn't the finite state machine familiar to programmers but a
mathematical model, coming from category theory, for a particularly
simple class of computing systems.  At a conceptual level, a finite
state machine can be considered a black box that receives a sequence
of input signal and produces an output signal for each input signal.
Internally, the black box maintains a *state* -- some sort of finite
summary representation of the sequence of input signals encountered so
far.  For each input signal, the box performs two operations.  In both
cases, the decision depends on the current internal state and the
identity of the input signal just received.

* **Chose next state**
* **Generate output token**

It is helpful to view the finite state machine layer as a mechanism to
simulate a *Markov chain* or *Markov process*.


Markov Chain for Discrete-Time Trajectories
=============================================

Roughly speaking, a *Markov chain*, :math:`\bf{X}`, is a probabilistic
system that makes random jumps among a finite set of distinct states,
:math:`s_0, s_1, s_2, \ldots, s_N` such that the probability of
choosing the next state, :math:`X_{n+1}` depends only on the current
state, :math:`X_n`.  In mathematical terms, the conditional
probabilities for state transitions must satisfy

.. math:: \mathcal{P}[X_{n+1} = s_{l} | X_0=s_i, X_1=s_j, \ldots, X_n=s_k] =
	  \mathcal{P}[X_{n+1} = s_{l} | X_{n}=s_k]

Since more distant history does not affect future behavior, Markov
chains are sometimes characterized as *memoryless*.

This relation can be iterated to compute
the conditional probabilities for multiple time steps

.. math:: \mathcal{P}[X_{n+2} = s_{m} | X_n=s_k] = \sum_{l} \mathcal{P}[X_{n+2} = s_{m} |
	  X_{n+1}=s_l] \mathcal{P}[X_{n+1} = s_{l} | X_{n}=s_k]

Note, the transition probabilities :math:`\mathcal{P}[X_{n+1} = s_{l} |
X_{n}=s_k]` may depend on time (through the index :math:`n`).  These so-called
time-inhomogeneous Markov chains arise when the system of interest is
driven by external entities.  Chains with time-independent conditional
transition probabilities are called time-homogeneous.  The dynamics of
a time-homogeneous Markov chain is completely determined by the
initial state and the transition probabilities.  All processes
considered in this document are time-homogeneous.

Markov Process for Continuous-Time Trajectories
=================================================

A *Markov process* is a generalization of
the Markov chain such that time is viewed as continuous rather than
discrete.  As a result, it makes sense to record the times at which
the transitions occur as part of the process itself.  

The first step in this generalization is to define a stochastic
process :math:`\bf{Y}` that includes the transition times as well as
the state, :math:`Y_{n} = (s_{j},t_{n})`.  

The second step is to treat time on a truly continuous basis by
defining a new stochastic process, :math:`\bf{Z}`, from :math:`\bf{Y}`
by the rule :math:`Z_{t} = s_k` in the time interval :math:`t_n \le t
< t_{n+1}` given :math:`Y_{n} = (s_k, t_n)` .  In other words,
:math:`\bf{Z}_{t}` is a piecewise constant version of :math:`\bf{Y}`
as shown in :ref:`piecewise_Z`

.. _piecewise_Z:

.. figure:: images/piecewise_Z.svg
   :scale: 100%
   :align: center

   Figure 4.  **Realization of a continuous time stochastic process and
   associated Markov chain.**

A realization of the process :math:`\bf{Y}` is defined by the closed
diamonds (left end points) alone.  Similarly, a realization of the
process :math:`\bf{Z}_t` is illustrated by the closed diamonds and
line segments.  The closed and open diamonds at the ends of the line
segment indicate that the segments include the left but not the right
end points.  

The memoryless property for Markov processes is considerably more
delicate than in the case of Markov chain because the time variable is
continuous rather than discrete.  In the case of :math:`\bf{Y}`, the
conditional probabilities for state transitions of must satisfy

.. math:: \mathcal{P}[Y_{n+1} = (s_{l},t_{n+1}) | Y_0=(s_i, t_0), Y_1=(s_j, t_1),
	  \ldots, Y_n=(s_k, t_n)] =
	  \mathcal{P}[Y_{n+1} = (s_{l}, t_{n+1}) | Y_{n}=(s_k, t_{n})]

The proper generalization of the requirement of time-homeogeneity
stated previously for Markov chains is that joint probability
be unchanged by uniform shifts in time

.. math:: \mathcal{P}[Z_{t+\tau} | Z_{s+\tau}] = \mathcal{P}[Z_{t} | Z_{s} ]

for :math:`0<s<t` and :math:`\tau > 0`.  Stochastic processes with
shift invariant state transition probabilities are called
*stationary*.  

When we examined hazard rates above, we were examining the rate
of transitions for a Markov process. The overall probability
of the next state of the Markov process is called the core
matrix,

.. math:: \mathcal{P}[Z_{t} | Z_{s} ]=Q_{ij}(t_{n+1}-t_n)

indicating a state change between the states :math:`(s_i,s_j)`.
The derivative of this is a rate,

.. math:: q_{ij}(t_{n+1}-t_n)=\frac{dQ_{ij}(t_{n+1}-t_n)}{dt},

which is a joint distribution over states and time intervals.
Normalization for this quantity sums over possible states
and future times,

.. math:: 1=\int_0^\infty \sum_j  q_{ij}(s)ds.

The survival, in terms of the core matrix, is

.. math::

   G_i(\tau)=1-\int_0^\tau \sum_k  q_{ik}(s)ds.

This means our hazard is

.. math::

   \lambda_{ij}(\tau)=\frac{q_{ij}(\tau)}{1-\int_0^\tau \sum_k  q_{ik}(s)ds}.

For the measles example, the set of future states :math:`j` of each individual
include symptomatic and all the possible other ways an individual
leaves the study, so you can think of :math:`j=\mbox{left town}`.
In practice, we build a hazard in two steps. First, count the probability
over all time for any one eventual state :math:`j`. This is the 
same stochastic probability :math:`\pi_{ij}` that is seen in Markov
chains. Second, measure the distribution of times at which 
intervals enter each new state :math:`j`, given that they are headed
to that state. This is called the holding time, :math:`h_{ij}(\tau)`,
and is a conditional probability. Together, these two give us
the core matrix,

.. math:: q_{ij}(\tau)=\pi_{ij}h_{ij}(\tau).

Note that :math:`h_{ij}(\tau)` is a density whose integral
:math:`H_{ij}(\tau)` is a cumulative distribution. If we write the
same equation in terms of probabilities, we see that it amounts
to separating the Markov process into a marginal and conditional
distribution.

.. math::

   \begin{eqnarray}
   q_{ij}(\tau)&=&\frac{d}{d\tau}P[Z_t|Z_s]\\
   &=&\frac{d}{d\tau}P[s_j|s_i,t_n]P[t_{n-1}-t_n\le\tau|s_i,s_j,t_n]\\
     & = & P[s_j|s_i,t_n]\frac{d}{d\tau}P[t_{n-1}-t_n\le\tau|s_i,s_j,t_n] \\
     & = & \pi_{ij}\frac{d}{d\tau}H_{ij}(\tau) \\
     & = & \pi_{ij}h_{ij}(\tau)
   \end{eqnarray}

Choosing the other option for the marginal gives us the
waiting time formulation for the core matrix. It corresponds
to asking first what is the distribution of times at which the
next event happens, no matter which event, and then asking
which events are more likely given the time of the event.

.. math::

   \begin{eqnarray}
   q_{ij}(\tau)&=&\frac{d}{d\tau}P[Z_t|Z_s]\\
   &=&\frac{d}{d\tau}P[s_j|s_i,t_n,t_{n+1}]P[t_{n-1}-t_n\le\tau|s_i,t_n]\\
     & = & \frac{d}{d\tau}(\Pi_{ij}(\tau)W_i(\tau)) \\
     & = & \pi_{ij}(\tau)\frac{d}{d\tau}W_i(\tau) \\
     & = & \pi_{ij}(\tau)w_{i}(\tau)
   \end{eqnarray}

While the waiting time density :math:`w_i(\tau)`, is the derivative
of the waiting time, we won't end up needing to relation
:math:`\pi_{ij}(\tau)` to :math:`\Pi_{ij}(\tau)` when finding trajectories
or computing hazards, so the more complicated relationship won't
be a problem.




Generalized Stochastic Petri Net for Bookkeeping
=================================================

A **generalized stochastic Petri net** (GSPN) is a formal way to
specify a a system of interacting, competing processes. Different
organisms can compete, but, for this system, the likelihood of
infecting a neighbor versus the likelihood of recovery are seen as
competing, as well.

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
