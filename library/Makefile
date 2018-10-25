# call the dependency generating tools to generate the dependency files,
# then call Makefile.1 to compile the libaray.

include ../Make.global_options

OBJ_EXT = .o
DIST_COMMON = Makefile Makefile.1 mkdep.pl mkfiles.pl
DIST_SOURCE = $(wildcard src/*.cpp)
DIST_HEADER = $(wildcard include/*.h)
DISTFILES = $(DIST_COMMON) $(DIST_SOURCE) $(DIST_HEADER)
mkinstalldirs = $(SHELL) ../mkinstalldirs

default: all

all:
	$(PERL) mkfiles.pl
	$(PERL) mkdep.pl
	cat make.ofiles | $(SED) -e "s/src\/\([^/]*\)\.o/obj\/\1$(OBJ_EXT)/g" > tmp.file
	mv tmp.file make.ofiles
	cat make.dep | $(SED) -e "s/^src/obj/g" > tmp.file
	mv tmp.file make.dep
	cat make.dep | $(SED) -e "s/\.o:/$(OBJ_EXT):/g" > tmp.file
	mv tmp.file make.dep
	mkdir -p obj
	mkdir -p lib
	$(MAKE) -f Makefile.1 default
	cd mpi && make

install:
	mkdir -p $(prefix)/include/AFEPack $(prefix)/include/AFEPack/mpi $(prefix)/lib
	cp -r $(AFEPACK_PATH)/library/include/*.h $(prefix)/include/AFEPack
	cd $(AFEPACK_PATH)/library/lib && cp *.so $(prefix)/lib
	cd mpi && make install
	echo "Please add the path $(prefix)/lib in LD_LIBRARY_PATH"
	echo "if the library can not be found automatically."

clean:
	-rm make.dep
	-rm make.files
	-rm make.ofiles
	-rm include/*~
	-rm src/*~
	-rm -rf obj
	-rm -rf lib
	cd mpi && make clean

distclean:

distdir: $(DISTFILES)
	@list='$(DISTFILES)'; for file in $$list; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  dir=`echo "$$file" | sed -e 's,/[^/]*$$,,'`; \
	  if test "$$dir" != "$$file" && test "$$dir" != "."; then \
	    dir="/$$dir"; \
	    $(mkinstalldirs) "$(distdir)$$dir"; \
	  else \
	    dir=''; \
	  fi; \
	  if test -d $$d/$$file; then \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -pR $(srcdir)/$$file $(distdir)$$dir || exit 1; \
	    fi; \
	    cp -pR $$d/$$file $(distdir)$$dir || exit 1; \
	  else \
	    test -f $(distdir)/$$file \
	    || cp -p $$d/$$file $(distdir)/$$file \
	    || exit 1; \
	  fi; \
	done


.PHONY: default all install clean distclean

