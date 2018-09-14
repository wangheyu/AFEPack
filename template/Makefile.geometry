# by R.Lie, Nov 1, 2002

include ../../Make.global_options

geometry = $(notdir $(PWD))
subdirs = $(wildcard *.[0-9]) $(wildcard *.[0-9].D)
csource = $(wildcard *.c)
cppsource = $(wildcard *.cpp)
source = $(csource) $(cppsource)
object = $(patsubst %.c, %.o, $(csource)) $(patsubst %.cpp, %.o, $(cppsource))
library = $(geometry).geometry.so

DIST_DATA = $(wildcard *.tmp_geo) $(wildcard *.crd_trs) $(wildcard *.out_nrm)
DISTFILES = $(DIST_DATA) $(source)
mkinstalldirs = $(SHELL) ../../../mkinstalldirs

all : default

default :
	$(MAKE) $(library)
	@$(foreach dir, $(subdirs), cd $(PWD)/$(dir) && $(RM) -f Makefile && $(LN_S) ../../Makefile.bas_fun ./Makefile && $(MAKE) ;)
	$(MAKE) symbolic

$(library) : $(object)
	$(CC) -shared -o $@ $(object)

symbolic :
	@$(foreach dir, $(subdirs), $(RM) -f $(dir).bas_fun ;)
	@$(foreach dir, $(subdirs), $(RM) -f $(dir).tmp_dof ;)
	@$(foreach dir, $(subdirs), $(RM) -f $(dir).bas_fun.so ;)
	@$(foreach dir, $(subdirs), $(LN_S) $(dir)/$(dir).bas_fun . ;)
	@$(foreach dir, $(subdirs), $(LN_S) $(dir)/$(dir).tmp_dof . ;)
	@$(foreach dir, $(subdirs), $(LN_S) $(dir)/$(dir).bas_fun.so . ;)

%.o : %.c
	$(CC) -c -o $@ $< $(CFLAGS) -fPIC

clean :
	-rm -f $(object)
	-rm -f $(library)
	@$(foreach dir, $(subdirs), $(RM) -f $(dir).bas_fun ;)
	@$(foreach dir, $(subdirs), $(RM) -f $(dir).tmp_dof ;)
	@$(foreach dir, $(subdirs), $(RM) -f $(dir).bas_fun.so ;)
	@$(foreach dir, $(subdirs), cd $(PWD)/$(dir) && $(MAKE) clean ;)
	-rm -f Makefile
	-rm -f *~

distdir:
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
	list='$(subdirs)'; for subdir in $$list; do \
	  if test "$$subdir" = .; then :; else \
	    test -d $(distdir)/$$subdir \
	    || mkdir $(distdir)/$$subdir \
	    || exit 1; \
	    (cd $$subdir && \
	      $(RM) -f Makefile && \
	      $(LN_S) ../../Makefile.bas_fun ./Makefile && \
	      $(MAKE) $(AM_MAKEFLAGS) \
	        top_distdir="$(top_distdir)" \
	        distdir=../$(distdir)/$$subdir \
	        distdir) \
	      || exit 1; \
	  fi; \
	done

.PHONY : default symbolic clean distdir

