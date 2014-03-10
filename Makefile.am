## Process this file with automake to produce Makefile.in
SEMIMARKOV_LIBRARY_VERSION=0:1:0
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}
AUTOMAKE_OPTIONS = foreign

AM_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB)

AM_CPPFLAGS = -std=c++11 -Isrc -Isrc/continuous -Isrc/gspn -Isrc/utility

#lib_LTLIBRARIES = semimarkov-0.1
#semimarkov_0_1_includedir = $(includedir)/semimarkov-0.1
#semimarkov_0_1_SOURCES =
#semimarkov_0_1_include_HEADERS = src/stochnet.h src/continuous/continuous_dynamics.h src/continuous/continuous_state.h src/continuous/marking.h src/gspn/distributions.h src/gspn/explicit_transitions.h src/gspn/petri_graph.h src/gspn/partial_core_matrix.h src/gspn/gspn.h src/utility/gspn_random.h src/utility/logging.h
#semimarkov_0_1_LDFLAGS = -version-info $(SEMIMARKOV_LIBRARY_VERSION)

bin_PROGRAMS = bvd sirmixed brownion

bvd_SOURCES = examples/bvd.cpp examples/cow_token.h src/stochnet.h src/continuous/continuous_dynamics.h src/continuous/continuous_state.h src/continuous/marking.h src/gspn/distributions.h src/gspn/explicit_transitions.h src/gspn/petri_graph.h src/gspn/partial_core_matrix.h src/gspn/gspn.h src/utility/gspn_random.h src/utility/logging.h

sirmixed_SOURCES = examples/sir_mixed.cpp src/stochnet.h src/continuous/continuous_dynamics.h src/continuous/continuous_state.h src/continuous/marking.h src/gspn/distributions.h src/gspn/explicit_transitions.h src/gspn/petri_graph.h src/gspn/partial_core_matrix.h src/gspn/gspn.h src/utility/gspn_random.h src/utility/logging.h

brownion_SOURCES = examples/brownion.cpp src/stochnet.h src/continuous/continuous_dynamics.h src/continuous/continuous_state.h src/continuous/marking.h src/gspn/distributions.h src/gspn/explicit_transitions.h src/gspn/partial_core_matrix.h src/gspn/gspn.h src/utility/gspn_random.h src/utility/logging.h

dist_noinst_SCRIPTS = autogen.sh
