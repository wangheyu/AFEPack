SUBDIRS = tetrahedron.0\
	tetrahedron.1\
	tetrahedron.1.D\
	tetrahedron.2\
	tetrahedron.2.D\
	tetrahedron.3.D\
	tetrahedron.edge.1


templatedir = $(includedir)/AFEPack/template/tetrahedron
template_LTLIBRARIES = libtetrahedron.geometry.la
libtetrahedron_geometry_la_SOURCES = tetrahedron.volume.c\
	tetrahedron.out_nrm.c\
	tetrahedron.crd_trs.c

EXTRA_DIST = tetrahedron.crd_trs tetrahedron.tmp_geo tetrahedron.out_nrm

install-data-local:
	cp $(srcdir)/tetrahedron.crd_trs\
	   $(srcdir)/tetrahedron.tmp_geo\
	   $(srcdir)/tetrahedron.out_nrm\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/tetrahedron.0/tetrahedron.0.bas_fun\
	   $(srcdir)/tetrahedron.0/tetrahedron.0.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/tetrahedron.1/tetrahedron.1.bas_fun\
	   $(srcdir)/tetrahedron.1/tetrahedron.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/tetrahedron.1.D/tetrahedron.1.D.bas_fun\
	   $(srcdir)/tetrahedron.1.D/tetrahedron.1.D.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/tetrahedron.2/tetrahedron.2.bas_fun\
	   $(srcdir)/tetrahedron.2/tetrahedron.2.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/tetrahedron.2.D/tetrahedron.2.D.bas_fun\
	   $(srcdir)/tetrahedron.2.D/tetrahedron.2.D.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/tetrahedron.3.D/tetrahedron.3.D.bas_fun\
	   $(srcdir)/tetrahedron.3.D/tetrahedron.3.D.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/tetrahedron.edge.1/tetrahedron.edge.1.bas_fun\
	   $(srcdir)/tetrahedron.edge.1/tetrahedron.edge.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/

uninstall-local:
	rm -f $(DESTDIR)$(templatedir)/tetrahedron.crd_trs\
              $(DESTDIR)$(templatedir)/tetrahedron.out_nrm\
	      $(DESTDIR)$(templatedir)/tetrahedron.tmp_geo
	rm -f $(DESTDIR)$(templatedir)/tetrahedron.0.bas_fun\
	   $(DESTDIR)$(templatedir)/tetrahedron.0.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/tetrahedron.1.bas_fun\
	   $(DESTDIR)$(templatedir)/tetrahedron.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/tetrahedron.1.D.bas_fun\
	   $(DESTDIR)$(templatedir)/tetrahedron.1.D.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/tetrahedron.2.bas_fun\
	   $(DESTDIR)$(templatedir)/tetrahedron.2.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/tetrahedron.2.D.bas_fun\
	   $(DESTDIR)$(templatedir)/tetrahedron.2.D.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/tetrahedron.3.D.bas_fun\
	   $(DESTDIR)$(templatedir)/tetrahedron.3.D.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/tetrahedron.edge.1.bas_fun\
	   $(DESTDIR)$(templatedir)/tetrahedron.edge.1.tmp_dof
