********************************
Requirements and Installation
********************************

This is a header-only library. Running ``./configure`` and ``make``
will build example codes and copy header files into the specified
``include`` directory.

Requirements
----------------

* C++ supporting the C++11 standard, which includes gcc 4.8 or newer
  and clang 3.4 or newer.

* Boost libraries, version 1.54 or newer.

* Make. The configure script and Makefile were created with autoconf
  and automake.

* Testing is regularly done on Linux and Mac, but the library is just C++
  so it should compile anywhere.


Installation
---------------

#. Determine the flavor of your boost libraries. The names of Boost
   binaries may contain extra information, most commonly ending in ``-mt``,
   such as libboost_system-mt.so. Locate your libraries with
   ``locate libboost_chrono`` or, on the mac, ``mdfind -name libboost_chrono``
   and ensure you have a full set of one of the
   variants. If your library sits somewhere like
   ``/usr/local/boost_1_54_0mt/lib/libboost_chrono-mt.so,`` then build
   with the following options to configure::

     ./configure --with-boost=/usr/local/boost_1_54_0mt --enable-variant=-mt

   Proper choice of boost directory and variant will lead to a string of
   ``yes`` in the configure output::

      checking for exit in -lpthread... yes
      checking for exit in -lboost_system... yes
      checking for exit in -lboost_random... yes
      checking for exit in -lboost_program_options... yes
      checking for exit in -lboost_filesystem... yes
      checking for exit in -lboost_date_time... yes
      checking for exit in -lboost_log... yes
      checking for exit in -lboost_log_setup... yes
      checking for exit in -lboost_unit_test_framework... yes


#. Modify the destination directory by adding a flag to configure::

     ./configure --prefix=$HOME

   With this option, the headers will install into the directory
   ``$HOME/include/semimarkov-0.1``.
   It is possible to modify the compiler and linker flags::

     ./configure CXX=g++ CXXFLAGS="-O2 -march=native"

   All of these sample flags can be used in combination, and help
   is found with ./configure --help.

#. Lastly, make the examples and install the headers::

      make
      make install

Should it be necessary to rebuild the ``configure`` script:

#. Run ``autoreconf --force --install``. This will run ``aclocal``,
   ``autoconf``, and ``automake``.