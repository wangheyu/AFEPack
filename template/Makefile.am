## Process this file with automake to create Makefile.in

DIST_COMMON = Makefile Makefile.in Makefile.am Makefile.bas_fun Makefile.geometry

all-am: symbolic Makefile

symbolic: 
	@echo "Make symbolic link Makefile for subdirectories." ; \
	list='$(SUBDIRS)'; \
	for subdir in $$list; do \
	  if test "$$subdir" != "."; then \
	    (cd $$subdir && $(LN_S) -f ../Makefile.geometry ./Makefile) ; \
	  fi; \
	done;

install:
	$(mkinstalldirs) $(prefix)/share/AFEPack
	cp -r ./* $(prefix)/share/AFEPack
	@echo "Remember to set the environment varible"
	@echo ""
	@echo "  AFEPACK_TEMPLATE_PATH"
	@echo ""
	@echo "when running your program. See document for details." 

SUBDIRS = . interval \
  triangle \
  twin_triangle \
  rectangle \
  tetrahedron \
  twin_tetrahedron \
  four_tetrahedron \
  quadrilateral \
  prism \
  hexahedron \
  recthexa \
  oblong \
  pyramid

