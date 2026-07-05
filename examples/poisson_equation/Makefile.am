exampledir = $(datadir)/@PACKAGE@/examples/poisson_equation

noinst_PROGRAMS = main
main_SOURCES = poisson_equation.cpp
main_LDADD = $(top_builddir)/src/libAFEPack.la

dist_example_DATA = poisson_equation.cpp poisson_equation.net README D.d
nodist_example_DATA = Makefile.sample run.sh.sample

EXTRA_DIST = Makefile.sample.in run.sh.sample.in
