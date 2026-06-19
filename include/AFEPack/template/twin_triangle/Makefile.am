SUBDIRS = twin_triangle.0\
	twin_triangle.1\
	twin_triangle.1.D\
	twin_triangle.2\
	twin_triangle.DG.1\
	twin_triangle.RT.1\
	twin_triangle.morley.1

templatedir = $(includedir)/AFEPack/template/twin_triangle
template_LTLIBRARIES = libtwin_triangle.geometry.la
libtwin_triangle_geometry_la_SOURCES = twin_triangle.volume.c\
	twin_triangle.crd_trs.c\
	twin_triangle.out_nrm.c\
	twin_triangle.to3d.crd_trs.c

EXTRA_DIST = twin_triangle.tmp_geo.new twin_triangle.iso.crd_trs twin_triangle.crd_trs twin_triangle.out_nrm twin_triangle.tmp_geo twin_triangle.to3d.crd_trs

install-data-local:
	cp $(srcdir)/twin_triangle.crd_trs\
	   $(srcdir)/twin_triangle.tmp_geo\
	   $(srcdir)/twin_triangle.out_nrm\
	   $(srcdir)/twin_triangle.to3d.crd_trs\
	   $(srcdir)/twin_triangle.tmp_geo.new\
	   $(srcdir)/twin_triangle.iso.crd_trs\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_triangle.0/twin_triangle.0.bas_fun\
	   $(srcdir)/twin_triangle.0/twin_triangle.0.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_triangle.1/twin_triangle.1.bas_fun\
	   $(srcdir)/twin_triangle.1/twin_triangle.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_triangle.1.D/twin_triangle.1.D.bas_fun\
	   $(srcdir)/twin_triangle.1.D/twin_triangle.1.D.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_triangle.2/twin_triangle.2.bas_fun\
	   $(srcdir)/twin_triangle.2/twin_triangle.2.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_triangle.DG.1/twin_triangle.DG.1.bas_fun\
	   $(srcdir)/twin_triangle.DG.1/twin_triangle.DG.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_triangle.morley.1/twin_triangle.morley.1.bas_fun\
	   $(srcdir)/twin_triangle.morley.1/twin_triangle.morley.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/twin_triangle.RT.1/twin_triangle.RT.1.bas_fun\
	   $(srcdir)/twin_triangle.RT.1/twin_triangle.RT.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/


uninstall-local:
	rm -f $(DESTDIR)$(templatedir)/twin_triangle.crd_trs\
              $(DESTDIR)$(templatedir)/twin_triangle.to3d.crd_trs\
              $(DESTDIR)$(templatedir)/twin_triangle.out_nrm\
	      $(DESTDIR)$(templatedir)/twin_triangle.tmp_geo\
	      $(DESTDIR)$(templatedir)/twin_triangle.tmp_geo.new\
	      $(DESTDIR)$(templatedir)/twin_triangle.iso.crd_trs
	rm -f $(DESTDIR)$(templatedir)/twin_triangle.0.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_triangle.0.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/twin_triangle.1.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_triangle.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/twin_triangle.1.D.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_triangle.1.D.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/twin_triangle.2.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_triangle.2.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/twin_triangle.DG.1.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_triangle.DG.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/twin_triangle.morley.1.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_triangle.morley.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/twin_triangle.RT.1.bas_fun\
	   $(DESTDIR)$(templatedir)/twin_triangle.RT.1.tmp_dof
