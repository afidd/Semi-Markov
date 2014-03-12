==============================
Requirements and Installation
==============================


Requirements:

* C++ supporting the C++11 standard, which includes gcc 4.8 or newer
  and clang 3.4 or newer.

* Boost libraries, version 1.54 or newer.

* Make. The configure script and Makefile were created with autoconnf
  and automake.


Installation:

#. Run "./configure --prefix=$HOME". Change the compiler or flags using
   environmental variables, CC, CFLAGS, and LDFLAGS before
   running configure.

#. Run "make". This makes examples only.

#. Run "make install". This will just copy header files to
   $HOME/includes.
