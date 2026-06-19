SUBDIRS = four_tetrahedron.0\
	four_tetrahedron.1\
	four_tetrahedron.1.D

templatedir = $(includedir)/AFEPack/template/four_tetrahedron
template_LTLIBRARIES = libfour_tetrahedron.geometry.la
libfour_tetrahedron_geometry_la_SOURCES = four_tetrahedron.volume.c\
	four_tetrahedron.out_nrm.c\
	four_tetrahedron.crd_trs.c

EXTRA_DIST = four_tetrahedron.crd_trs four_tetrahedron.tmp_geo four_tetrahedron.out_nrm

install-data-local:
	cp $(srcdir)/four_tetrahedron.crd_trs\
	   $(srcdir)/four_tetrahedron.tmp_geo\
	   $(srcdir)/four_tetrahedron.out_nrm\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/four_tetrahedron.0/four_tetrahedron.0.bas_fun\
	   $(srcdir)/four_tetrahedron.0/four_tetrahedron.0.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/four_tetrahedron.1/four_tetrahedron.1.bas_fun\
	   $(srcdir)/four_tetrahedron.1/four_tetrahedron.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/four_tetrahedron.1.D/four_tetrahedron.1.D.bas_fun\
	   $(srcdir)/four_tetrahedron.1.D/four_tetrahedron.1.D.tmp_dof\
	   $(DESTDIR)$(templatedir)/

uninstall-local:
	rm -f $(DESTDIR)$(templatedir)/four_tetrahedron.crd_trs\
              $(DESTDIR)$(templatedir)/four_tetrahedron.out_nrm\
	      $(DESTDIR)$(templatedir)/four_tetrahedron.tmp_geo
	rm -f $(DESTDIR)$(templatedir)/four_tetrahedron.0.bas_fun\
	   $(DESTDIR)$(templatedir)/four_tetrahedron.0.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/four_tetrahedron.1.bas_fun\
	   $(DESTDIR)$(templatedir)/four_tetrahedron.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/four_tetrahedron.1.D.bas_fun\
	   $(DESTDIR)$(templatedir)/four_tetrahedron.1.D.tmp_dof
