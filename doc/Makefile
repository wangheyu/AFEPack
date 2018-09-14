# by R.Lie, Nov 1, 2002

include ../Make.global_options

default : all

ifeq (x$(DOCGEN), x)
  all :
	echo "doxygen is not found that not able to generate document"
else
  all : $(AFEPACK_PATH)/library/include/*.h
	$(DOCGEN) AFEPack.cfg
endif

FIGURES = $(wildcard graphics/*jpg)
DISTFILES = Makefile AFEPack.cfg \
  index.html index-1.html index-2.html \
  index-3.html index-4.html index-5.html \
  index-6.html install.texi $(FIGURES)
  
mkinstalldirs = $(SHELL) ../mkinstalldirs

install:
	install-dir='$(prefix)/doc/$(PACKAGE)-$(VERSION)'; \
	mkdir $(install-dir) ; \
	mkdir $(install-dir)/graphics; \
	@list='$(DISTFILES)'; for file in $$list; do \
	  cp "$$file" $(install-dir); \
	done; \
	cp -R html $(install-dir)/share/doc/AFEPack
	cp -R man $(install-dir)/man

clean:
	-rm -rf html
	-rm -rf latex
	-rm -rf man
	-rm -rf rtf

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

.PHONY : default all clean distdir
