SUBDIRS = twin_tetrahedron.0\
	twin_tetrahedron.1\
	twin_tetrahedron.1.D

templatedir = $(includedir)/AFEPack/template/twin_tetrahedron
template_LTLIBRARIES = libtwin_tetrahedron.geometry.la
libtwin_tetrahedron_geometry_la_SOURCES = twin_tetrahedron.volume.c\
	twin_tetrahedron.out_nrm.c\
	twin_tetrahedron.crd_trs.c

EXTRA_DIST = twin_tetrahedron.crd_trs twin_tetrahedron.tmp_geo twin_tetrahedron.out_nrm


install-data-local:
	cp $(srcdir)/twin_tetrahedron.crd_trs\
	   $(srcdir)/twin_tetrahedron.tmp_geo\
	   $(srcdir)/twin_tetrahedron.out_nrm\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_tetrahedron.0/twin_tetrahedron.0.bas_fun\
	   $(srcdir)/twin_tetrahedron.0/twin_tetrahedron.0.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_tetrahedron.1/twin_tetrahedron.1.bas_fun\
	   $(srcdir)/twin_tetrahedron.1/twin_tetrahedron.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_tetrahedron.1.D/twin_tetrahedron.1.D.bas_fun\
	   $(srcdir)/twin_tetrahedron.1.D/twin_tetrahedron.1.D.tmp_dof\
	   $(DESTDIR)$(templatedir)/

uninstall-local:
	rm -f $(DESTDIR)$(templatedir)/twin_tetrahedron.crd_trs\
              $(DESTDIR)$(templatedir)/twin_tetrahedron.out_nrm\
	      $(DESTDIR)$(templatedir)/twin_tetrahedron.tmp_geo
	rm -f $(DESTDIR)$(templatedir)/twin_tetrahedron.0.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_tetrahedron.0.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/twin_tetrahedron.1.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_tetrahedron.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/twin_tetrahedron.1.D.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_tetrahedron.1.D.tmp_dof
