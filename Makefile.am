## Process this file with automake to produce Makefile.in
SEMIMARKOV_LIBRARY_VERSION=0:1:0
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}
AUTOMAKE_OPTIONS = foreign

AM_LDFLAGS = $(BOOST_LDFLAGS) -fPIC

AM_CPPFLAGS = -std=c++11 -Isrc -Isrc/continuous -Isrc/gspn -Isrc/utility $(BOOST_CPPFLAGS) -fPIC

bin_PROGRAMS = weiss bvd sirmixed meta disttest sirexp

semimarkov_0_1_includedir = $(includedir)/semimarkov-0.1
semimarkov_0_1_include_HEADERS = src/stochnet.hpp src/smv.hpp src/smv_config.hpp src/continuous/continuous_dynamics.hpp src/continuous/continuous_state.hpp src/continuous/local_marking.hpp src/continuous/marking.hpp src/continuous/stochastic_dynamics.hpp src/gspn/build_graph.hpp src/gspn/distributions.hpp src/gspn/explicit_transitions.hpp src/gspn/finite_state_machine.hpp src/gspn/gspn.hpp src/gspn/partial_core_matrix.hpp src/gspn/petri_graph.hpp src/utility/gspn_random.hpp src/utility/logging.hpp  src/utility/smv_algorithm.hpp
#semimarkov_0_1_LDFLAGS = -version-info $(SEMIMARKOV_LIBRARY_VERSION)


bvd_SOURCES = example/bvd.cpp example/cow_token.hpp $(semimarkov_0_1_include_HEADERS)

sirmixed_SOURCES = example/sir_mixed.cpp $(semimarkov_0_1_include_HEADERS)

sirexp_SOURCES = example/sir_exp.cpp $(semimarkov_0_1_include_HEADERS)

meta_SOURCES = example/meta.cpp $(semimarkov_0_1_include_HEADERS)

weiss_SOURCES = example/Weiss/weiss.cpp example/Weiss/weiss.hpp $(semimarkov_0_1_include_HEADERS)

disttest_SOURCES = test/distributions_test.cpp src/stochnet.hpp src/gspn/distributions.hpp src/utility/gspn_random.hpp src/utility/logging.hpp

dist_noinst_SCRIPTS = autogen.sh
