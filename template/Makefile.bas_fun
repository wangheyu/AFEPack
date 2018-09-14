# by R.Lie, Nov 1, 2002

include ../../../Make.global_options

csource = $(wildcard *.c)
cppsource = $(wildcard *.cpp)
object = $(patsubst %.c, %.o, $(csource)) \
         $(patsubst %.cpp, %.o, $(cppsource))
library = $(patsubst %.o, %.so, $(object))

DIST_DATA = $(wildcard *.bas_fun) $(wildcard *.tmp_dof)
DISTFILES = $(DIST_DATA) $(csource) $(cppsource)
mkinstalldirs = $(SHELL) ../../../../mkinstalldirs

%.o : %.c
	$(CC) -c $< -I$(AFEPACK_PATH)/library/include $(CFLAGS)

%.o : %.cpp
	$(CXX) -c $< -I$(AFEPACK_PATH)/library/include $(CXXFLAGS)

%.so : $(object)
	$(CC) -shared -o $@ $(object)

all : default

default : $(library)

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

clean :
	-rm -f $(object)
	-rm -f $(library)
	-rm -f Makefile
	-rm -f *~

.PHONY : default clean distdir

