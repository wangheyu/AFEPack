exampledir = $(datadir)/doc/AFEPack/example/DysonSphere

dist_example_DATA = main.cpp Makefile.sample run.sh.sample showMesh.net README Cube.geo

install-data-hook:
	chmod 666 $(DESTDIR)$(exampledir)/Makefile.sample
	sed -i '/CPPFLAGS =/ c CPPFLAGS=-I$(prefix)/include' $(DESTDIR)$(exampledir)/Makefile.sample
	sed -i '/LDFLAGS =/ c LDFLAGS=-L$(prefix)/lib' $(DESTDIR)$(exampledir)/Makefile.sample
	chmod 777 $(DESTDIR)$(exampledir)/run.sh.sample
	sed -i '/export LD_LIBRARY_PATH/ c export LD_LIBRARY_PATH=$(prefix)/lib' $(DESTDIR)$(exampledir)/run.sh.sample
	sed -i '/export AFEPACK_PATH/ c export AFEPACK_PATH=$(prefix)/include/AFEPack' $(DESTDIR)$(exampledir)/run.sh.sample
