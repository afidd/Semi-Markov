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
for knowledge that it has not yet fired. This renormalization
is apparently less common, because it is not in the API of
Boost::Random or std::random.

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

**ExponentialDistribution** This is the classic Markovian distribution.

**BoostDistribution** This takes any Boost distribution and converts it
to account for the renormalization of a time shift.

.. math::

   quantile(distribution, 1-(1-rand(rng))/normal)
   *(1-cdf(distribution, t_c-t_e))-(t_c-t_e)

This transformation accounts for the time shifting.

