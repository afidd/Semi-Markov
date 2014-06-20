################################
Semi-Markov Library FAQ
################################



Testing
^^^^^^^^^^^

There are three levels of testing.

* **Build Tests** Regular testing is done on Linux and Mac with gcc 4.7, gcc 4.8, clang 3.3,
  and clang 3.4. There is a `Jenkins <http://jenkins-ci.org/>`_ instance locally that runs
  builds inside of `Docker <http://docker.io/>`_ containers. This handles builds on Ubuntu
  and Fedora. Mac builds are tested by hand.
  A description of how to create Docker instances, and the script for these,
  is the in ``test`` subdirectory.

.. figure:: images/jenkins_main.*
   :align: center
   :width: 400px

* **System-level Tests** There are currently just a few of these, namely a set of SIR
  models whose theoretical properties are known, contained in examples/sir_mixed.cpp.
  We verify the percentage of runs which become epidemic.

* **Unit Tests** The distributions (exponential, Weibull, and others) have a unit test
  called disttest, built from the test directory. This samples the distributions
  and compares Kolmogorov-Smirnov estimators of those samples with theoretical results.


Code Style Guide
^^^^^^^^^^^^^^^^^^^^

The code follows the `Boost library guidelines <http://www.boost.org/development/requirements.html>`_ for compatibility, but Boost isn't
a style guide, so it follows the
`Google C++ Style Guide <http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml>`_.


Repository
^^^^^^^^^^^^^^

The repository on Github uses the rebase workflow from
`A Rebase Workflow for Git <http://randyfay.com/content/rebase-workflow-git>`_.
The public branch is ``master`` and topic branches are named and rebased
as we work.


