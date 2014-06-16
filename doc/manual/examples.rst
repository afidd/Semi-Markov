*********************
Examples
*********************


Weiss Two-state Brownion
-------------------------
The `Weiss two state Brownion example`_ is in a separate PDF.

.. _Weiss two state Brownion example: weiss.pdf

Susceptible-infected-susceptible (SIS) model of infectious disease transmission
-------------------------------------------------------------------------------


Management of dairy herds
-------------------------

As an example, let's look at a simple model for the farm management of
dairy cows, following [Viet:2004]_.  There are four groups, calves,
heifers before breeding (heifer 1), heifers ready for breeding (heifer
2), and dairy cows which give birth to new calves. We might make a
sketch as shown here.

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
