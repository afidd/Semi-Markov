==============================
Requirements and Installation
==============================


Requirements:

* C++ supporting the C++11 standard, which includes gcc 4.8 or newer
  and clang 3.4 or newer.

* Boost libraries, version 1.54 or newer.

* Make. The configure script and Makefile were created with autoconnf
  and automake.

* The random number generation can use std::random, Boost::random,
  or RNGSSELIB or PRAND from the `Computer Physics Communications Program Library of ScienceDirect <http://cpc.cs.qub.ac.uk/overview.html>`_.
  RNGSSELIB uses SSE instructions and more modern generators. PRAND
  augments this capability with GPU generation of random numbers.

Installation:

#. Run "./configure --prefix=$HOME". Change the compiler or flags using
   environmental variables, CC, CFLAGS, and LDFLAGS before
   running configure. The only real option is specification of the 
   location of Boost and the variant to link. If Boost libraries are in
   /usr/local/lib/libboost and headers are in /usr/local/include/boost/,
   then specify "./configure --with-boost=/usr/local". Boost has library
   variants. If your libraries look like "libboost_system-mt.so", then
   use "./configure --enable-variant=-mt".

#. Run "make". This makes examples only.

#. Run "make install". This will just copy header files to
   $HOME/includes.
