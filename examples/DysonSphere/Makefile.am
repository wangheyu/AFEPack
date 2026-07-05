exampledir = $(datadir)/@PACKAGE@/examples/DysonSphere

noinst_PROGRAMS = main
main_SOURCES = main.cpp
main_LDADD = $(top_builddir)/src/libAFEPack.la

dist_example_DATA = main.cpp showMesh.net README Cube.geo
nodist_example_DATA = Makefile.sample run.sh.sample

EXTRA_DIST = Makefile.sample.in run.sh.sample.in

# install-data-hook:
# 	chmod 666 $(DESTDIR)$(exampledir)/Makefile.sample
# 	sed -i '/CPPFLAGS =/ c CPPFLAGS=-I$(prefix)/include' $(DESTDIR)$(exampledir)/Makefile.sample
# 	sed -i '/LDFLAGS =/ c LDFLAGS=-L$(prefix)/lib' $(DESTDIR)$(exampledir)/Makefile.sample
# 	chmod 777 $(DESTDIR)$(exampledir)/run.sh.sample
# 	sed -i '/export LD_LIBRARY_PATH/ c export LD_LIBRARY_PATH=$(prefix)/lib' $(DESTDIR)$(exampledir)/run.sh.sample
# 	sed -i '/export AFEPACK_PATH/ c export AFEPACK_PATH=$(prefix)/include/AFEPack' $(DESTDIR)$(exampledir)/run.sh.sample
