==================================================
Introduction
==================================================

This is a library to simulate trajectories of semi-Markov processes in
continuous time. It represents these processes as Generalized
Stochastic Petri Nets (GSPN). The goal of this library is to create a 
clear specification for, and efficient computation of, these systems.



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
