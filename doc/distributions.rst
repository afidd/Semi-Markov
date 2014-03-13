==========================
Continuous Distributions
==========================

The core matrix of a semi-Markov system is defined on a joint
discrete and continuous space, where the discrete set is of
next possible states of the system and the continuous space
is of the next times when the system arrives at that state.
A competing process description of a semi-Markov system,
describes decides which state is next and when that state
happens by asking which of a set of continuous distributions,
one per competing process, is next to fire.

The GSPN, in particular, decides that, at the moment a transition
first becomes enabled by the marking, it chooses the distribution
of firing times. At any later time when a different, independent,
transition has fired without affecting the first, the first's
distribution doesn't change. It just gets renormalized to account
for knowledge that it has not yet fired. If we label the 
cumulative distribution by :math:`F(t)` and a distribution is shifted
by a time :math:`\Delta=t_0-t_e`, then the new distribution is
sampled with a shifted quantile,

.. math::

  t-t_0=Q(U,\Delta)=F^{-1}(U(1-F(\Delta))+F(\Delta))-\Delta.

Given any distributions from Boost::Random or std::random, this
formula calculates the shifted quantile from the original quantile
and cumulative distribution.

Transitions for this library implement the following interface.::

   template<typename RNG>
   class TransitionDistribution
   {
   public:
       virtual double sample(double enabling_time, double current_time, RNG& rng);
   };

Each sample can be taken from the difference between the two times,
:math:`(t_c-t_e)`, but passing both separately will be better for vectorization
of sampling.

The following distributions are currently in the library.

**NoDistribution**
This has no parameters and always returns infinity.

**ExponentialDistribution** The `exponential distribution <http://en.wikipedia.org/wiki/Exponential_distribution>`_
is the classic Markovian distribution.

.. math::

  F(x) = 1-e^{-λx}

  Q(x) = -(1/\lambda)\ln(1-x)

  Q(x,\Delta)=Q(x)

**WeibullDistribution** `Weibull distributions <http://en.wikipedia.org/wiki/Weibull_distribution>`_ can model either
infant mortality or aging processes.

.. math::

  F(x)=1-e^{-\left(x/\lambda\right)^k}

  Q(p; k,\lambda)=\lambda\left[-\ln(1-p)\right]^{1/k}

  Q(p,\Delta; k,\lambda)=\lambda\left[-\ln(1-p)+\left(\Delta/\lambda\right)^k\right]^{1/k}-\Delta

This transformation accounts for the time shifting.

**GammaDistribution** This uses the Boost::Math::gamma_distribution.
It has two parameters, shape and scale.