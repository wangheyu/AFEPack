#exampledir = $(docdir)/example/poisson_equation

exampledir = $(datadir)/doc/AFEPack/example/local_refine/2d

dist_example_DATA = 2d.cpp  2d.net  Coarse.d  D.d  Makefile.sample  README  run.sh.sample

install-data-hook:
	chmod 666 $(DESTDIR)$(exampledir)/Makefile.sample
	sed -i '/CPPFLAGS =/ c CPPFLAGS=-I$(prefix)/include' $(DESTDIR)$(exampledir)/Makefile.sample
	sed -i '/LDFLAGS =/ c LDFLAGS=-L$(prefix)/include' $(DESTDIR)$(exampledir)/Makefile.sample
	chmod 777 $(DESTDIR)$(exampledir)/run.sh.sample
	sed -i '/export LD_LIBRARY_PATH/ c export LD_LIBRARY_PATH=$(prefix)/lib' $(DESTDIR)$(exampledir)/run.sh.sample
	sed -i '/export AFEPACK_PATH/ c export AFEPACK_PATH=$(prefix)/include/AFEPack' $(DESTDIR)$(exampledir)/run.sh.sample
