************************************
Continuous Distributions Reference
************************************

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
sampled with a shifted quantile. If we know the original cumulative
distribution function, :math:`F(t,t_e),` where :math:`t_e` is the
enabling time of the distribution function, then we can calulate
the distribution, *given that it hasn't fired before time* :math:`t_0`,
with

.. math::

  F(t, t_0, t_e)=\frac{F(t,t_e)-F(t_0,t_e)}{1-F(t_0,t_e)}

This makes the calculation of the quantile, using the original
quantile, :math:`F^{-1},` look like

.. math::

  t-t_0=Q(U,\Delta)=F^{-1}(U(1-F(\Delta))+F(\Delta))-\Delta.

Given any distributions from Boost::Random or std::random, this
formula calculates the shifted quantile from the original quantile
and cumulative distribution.

TransitionDistribution Base Class
====================================

Transitions for this library implement the following interface.

.. cpp:class:: afidd::smv::TransitionDistribution<RandomGenerator>

   This is the pure virtual base class for distributions.

.. cpp:function:: double Sample(double current_time, RandomGenerator& rng) const

   Return a sample from the random variable in absolute time
   after the given current time.

.. cpp:function:: double EnablingTime() const

   Each distribution remembers its enabling time in absolute time
   since the start of the simulation.
 
.. cpp:function:: bool BoundedHazard() const

   Some distributions have a hazard function, defined as the ratio
   of the probability density to the survival, which is bounded for
   finite times. These functions will have a well-defined hazard
   integral. This method tells the caller whether the hazard integral
   is a reasonable approach to sample the distribution.

.. cpp:function:: double HazardIntegral(double t0, double t1) const

   This returns the integral of the hazard between absolute time
   :math:`t_0` and :math:`t_1`.

.. cpp:function:: double ImplicitHazardIntegral(double xa, double t0) const

   This implicitly solves for a quantile by integrating the hazard.
   For certain distributions, this is analytic. The function returns
   :math:`t` in

.. math::
   
   x_a=\int_{t_0}^{t} \lambda(s) ds.

The template argument RandomGenerator fulfills the concept of the Boost
random number generator concept and the std::random random number
generator concept.
Each sample can be taken from the difference between the two times,
:math:`(t_c-t_e)`, but passing both separately will be better for vectorization
of sampling.

The following distributions are currently in the library.

Exponential Distribution
===========================

.. cpp:class:: afidd::smv::ExponentialDistribution<RandomNumberGenerator>

   The `exponential distribution <http://en.wikipedia.org/wiki/Exponential_distribution>`_
   is the classic Markovian distribution. The cumulative
   distribution and quantile are defined as the following.

.. math::

  F(x) = 1-e^{-λx}

  Q(x) = -(1/\lambda)\ln(1-x)

  Q(x,\Delta)=Q(x)


.. cpp:function:: afidd::smv::ExponentialDistribution::ExponentialDistribution(double lambda, double enabling_time)

   The constructor takes the parameter `lambda,` an enabling time for the
   distribution as an absolute system time.


.. cpp:class:: afidd::smv::ShiftedExponentialDistribution<RandomNumberGenerator>

   The `exponential distribution <http://en.wikipedia.org/wiki/Exponential_distribution>`_
   is the classic Markovian distribution.
   The shift is a displacement
   of the cumulative distribution function by an amount :math:`t_s.`

.. math::

     F(x) = 1-e^{-λ(x-t_s)}


.. cpp:function:: afidd::smv::ShiftedExponentialDistribution::ExponentialDistribution(double lambda, double enabling_time, double shift=0.0, double normal=1.0)

   The constructor takes the parameter `lambda,` an enabling time for the
   distribution as an absolute system time, and a normalization constant.
   If `normal` is less than one, it represents the probability that
   this distribution will fire at all. 


Weibull Distribution
=======================

.. cpp:class:: afidd::smv::WeibullDistribution<RandomNumberGenerator>

   `Weibull distributions <http://en.wikipedia.org/wiki/Weibull_distribution>`_ can model either infant mortality or aging processes. Parameters
   may be defined different ways. This class uses the following
   cumulative distribution and quantile.

.. math::

      F(x)=1-e^{-\left(x/\lambda\right)^k}

      Q(p; k,\lambda)=\lambda\left[-\ln(1-p)\right]^{1/k}

      Q(p,\Delta; k,\lambda)=\lambda\left[-\ln(1-p)+\left(\Delta/\lambda\right)^k\right]^{1/k}-\Delta

.. cpp:function:: afidd::smv::WeibullDistribution::WeibullDistribution(double lambda, double k, double enabling_time, double shift, double normal=1.0)

   This creates a Weibull distribution with parameters as defined above.
   The shift moves the distribution to the right.

Gamma Distribution
===========================

.. cpp:class:: afidd::smv::GammaDistribution<RandomNumberGenerator>

   This uses the Boost::Math::gamma_distribution.
   It has two parameters, shape and scale.

.. cpp:function:: afidd::smv::GammaDistribution::GammaDistribution(double alpha, double theta, double enabling_time, double shift=0.0, double normal=1.0)

   The constructor initializes the two parameters, :math:`\alpha` and :math:`\theta.` It also sets the enabling time and optional shift and normal.

Piecewise Linear Distribution
===============================

.. cpp:class:: afidd::smv::PiecewiseLinearDistribution<RandomNumberGenerator>

   This distribution represents piecewise, linear, continuous distributions.
   It is an expansion on the `std::piecewise_linear_distribution` from
   the `std::random` header. The piecewise curve defines an un-normalized
   probability density function, from which the cumulative distribution
   function is calculated.

.. cpp:function:: afidd::smv::PiecewiseLinearDistribution::PiecewiseLinearDistribution( const std::vector<double>& b, const std::vector<double>& w, double enabling_time, double shift=0.0, double normal=1.0)

   The vector `b` specifies intercepts on the x-axis. The domain of the
   probability distribution function is from the first to last value of
   `b`. The weight vector, `w,` is the height of the unnormalized
   function at each point `b.` The arrays `b` and `w` must have at
   least two points and must be the same length.
   The `shift` moves the whole distribution to the right. `normal`
   is the probability that this distribution will fire at all.
