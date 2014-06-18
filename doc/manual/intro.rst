*********************
Introduction
*********************

The ``Semi-Markov`` library is designed to streamline the process of
creating efficient simulations of a large class of systems called
semi-Markov processes [Howard:1971]_.  Semi-Markov processes naturally
arise in many contexts including epidemiology [Viet:2004]_, physiology,
ecology, atmospheric sciences, reliability engineering and risk
management.  This broad range of applications suggest the value of
designing a generic library for simulating complex semi-Markov
processes, independent of the particular application area.  

The unifying idea on which this library is based is that there
typically are many different pathways for a complex system to evolve
between timesteps.  Each pathway can be viewed as an elementary
stochastic process with a user specified time-dependent transition
rates and a rule for modifying the overall internal state of the
system.  At each instant of time, these elementary processes
"compete", figuratively speaking, for the chance to change the state
of the whole system.  Each time step in the simulation corresponds to
an event -- a "winner" is selected thus changing the internal state of
the system and the sampling from the corresponding statistical
distribution to determine the time increment.  This competing process
view provides a framework for users to develop simulations for complex
models in an incremental manner.

It is easy to show that competing processes with exponentially
distributed transition times have time-independent transition rates.
This is the norm in some application areas such as chemical kinetics.
In contrast, it is manifestly inappropriate for many biological
applications such as physiology, ecology and epidemiology.  For
example, a classic paper by Stocks [Stocks:1931]_ clearly shows that
the latent period for measles (the distribution times between
infection and the appearance of symptoms) does not follow an
exponential distribution. Stocks' raw data from cases in London circa
1931, along with optimal fits to exponential, gamma, Weibull, and
log-normal distributions computed using the ``SciPy`` statistical
library, are shown in :ref:`latent_period`.  The fit of the data to
the exponential distribution is very poor while the fits to the other
distributions are very good.

.. _latent_period:

.. figure:: images/Stocks_fitted_data.*
   :scale: 50%
   :align: center

   Figure 1.  Distribution of latent periods for measles in London
   circa 1931

This simple example shows that exclusive reliance on exponential
distributions may lead systematic biases in stochastic simulations of
epidemiological process.  Therefore, this library provides support for
general semi-Markov models based on competing processes with general
probability distributions of transition times.  


Organization of library
---------------------------

It is implemented using three cooperating layers:

* **Finite state machine**: High-level interface for initializing the
  system, iterating over time steps and gathering relevant tracing
  data for post-processing.

* **Semi-Markov process**: "Middleware" responsible for statistically
  unbiased choice among all possible competing processes at a given
  time step.

* **Generalized stochastic Petri net**: Low-level coordination and
  bookkeeping related to the user-defined competing processes
  including distributions of transition times, modification of system
  state and various dependence relationships.

.. figure:: images/large-structureb.*
   :align: center
   :width: 400px

This organization has many practical advantages:

* The semi-Markov process layer can be viewed as a very efficient,
  general purpose, stochastic simulation engine that supports
  arbitrary statistical distributions for event times.  This layer
  contains no model-specific user code, thus can be independently
  verified and validated.

* Typically, the itself model is completely defined by instantiating
  components of the Petri net.  The Petri net layer automates most of
  the tedious and error-prone bookkeeping steps associated with the
  execution of the model.

* The library strictly enforces a separation of the static components
  that define the structural aspects of the model and the dynamic
  components that define the evolving state during a simulation.  This
  separation makes it possible to detect many critical programming
  errors associated with multithreading at compile time.


Larger Context
-----------------
Making a dynamical model, such as this library creates, is one
step in achieving larger goals. The goal could be *model selection*
to ask which conceptual process best describes observed data.
It could be *optimization* of interventions to hinder or encourage
an outcome. For example, see [Hartig2011]_ on using dynamical models
for *inference.* These larger goals guide construction of a
conceptual model for the problem, one that includes those behaviors
that seem most relevant.

The steps from conception to having a running program which uses
the Semi-Markov library are:

.. figure:: images/multi-large.*
   :align: center
   :width: 800px

1. A conceptual model for the system. What actors and resources are
   in the domain? What behaviors are important? A paper by
   [Grimm2010]_ and others goes into more detail for the ecology realm.

2. A survival analysis of that conceptual model. There are books and
   statistics courses on survival analysis. Field observations are
   processed with statistical estimators to produce estimated
   hazards for transition.

3. Now write C++ code, using the Semi-Markov library API, to
   record the states of the system and hazards for transition
   among those states. The data structure holding this information
   is called the GSPN.

4. When the library compiles, it generates, from the GSPN provided,
   an event loop to yield the next state change in the system and
   the time at which that change occurs.

5. Most often dynamical models build with this library will be
   part of the inner loop of a larger mission, for instance
   a Bayesian MCMC or an optimization. In the least, an ensemble of
   trajectories generated by a dynamical model will be summarized
   with statistical estimators at the finish.


Acknowledgements
--------------------

This library was created by the Analytical Framework for Infectious
Disease Dynamics (AFIDD) group at Cornell University in conjunction
with the USDA Agricultural Research Service.  This work was supported by the Science \& Technology Directorate, Department of Homeland Security via interagency agreement no. HSHQDC-10-X-00138.


Availability and distribution
-------------------------------

This library is in the public domain.  


