templatedir = $(includedir)/AFEPack/template/triangle

template_LTLIBRARIES = libtriangle.geometry.la

libtriangle_geometry_la_SOURCES = triangle.volume.c\
	triangle.crd_trs.c\
	triangle.to3d.crd_trs.c\
	triangle.out_nrm.c

SUBDIRS = triangle.0\
	triangle.1\
	triangle.1.D\
	triangle.2\
	triangle.3\
	triangle.4\
	triangle.DG.1\
	triangle.DG.2\
	triangle.DG.3\
	triangle.RT.1\
	triangle.RT.2\
	triangle.dual.0\
	triangle.edge.1\
	triangle.edge.Full.1\
	triangle.morley.1

EXTRA_DIST = triangle.crd_trs\
	triangle.to3d.crd_trs\
	triangle.out_nrm\
	triangle.tmp_geo

install-data-local:
	cp $(srcdir)/triangle.crd_trs\
	   $(srcdir)/triangle.tmp_geo\
	   $(srcdir)/triangle.out_nrm\
	   $(srcdir)/triangle.to3d.crd_trs\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.0/triangle.0.bas_fun\
	   $(srcdir)/triangle.0/triangle.0.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.1/triangle.1.bas_fun\
	   $(srcdir)/triangle.1/triangle.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.1.D/triangle.1.D.bas_fun\
	   $(srcdir)/triangle.1.D/triangle.1.D.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.2/triangle.2.bas_fun\
	   $(srcdir)/triangle.2/triangle.2.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.3/triangle.3.bas_fun\
	   $(srcdir)/triangle.3/triangle.3.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.4/triangle.4.bas_fun\
	   $(srcdir)/triangle.4/triangle.4.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.DG.1/triangle.DG.1.bas_fun\
	   $(srcdir)/triangle.DG.1/triangle.DG.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.DG.2/triangle.DG.2.bas_fun\
	   $(srcdir)/triangle.DG.2/triangle.DG.2.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.DG.3/triangle.DG.3.bas_fun\
	   $(srcdir)/triangle.DG.3/triangle.DG.3.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.dual.0/triangle.dual.0.bas_fun\
	   $(srcdir)/triangle.dual.0/triangle.dual.0.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.edge.1/triangle.edge.1.bas_fun\
	   $(srcdir)/triangle.edge.1/triangle.edge.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.edge.Full.1/triangle.edge.Full.1.bas_fun\
	   $(srcdir)/triangle.edge.Full.1/triangle.edge.Full.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.morley.1/triangle.morley.1.bas_fun\
	   $(srcdir)/triangle.morley.1/triangle.morley.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.RT.1/triangle.RT.1.bas_fun\
	   $(srcdir)/triangle.RT.1/triangle.RT.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/triangle.RT.2/triangle.RT.2.bas_fun\
	   $(srcdir)/triangle.RT.2/triangle.RT.2.tmp_dof\
	   $(DESTDIR)$(templatedir)/

uninstall-local:
	rm -f $(DESTDIR)$(templatedir)/triangle.crd_trs\
              $(DESTDIR)$(templatedir)/triangle.to3d.crd_trs\
              $(DESTDIR)$(templatedir)/triangle.out_nrm\
	      $(DESTDIR)$(templatedir)/triangle.tmp_geo
	rm -f $(DESTDIR)$(templatedir)/triangle.0.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.0.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.1.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.1.D.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.1.D.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.2.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.2.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.3.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.3.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.4.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.4.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.DG.1.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.DG.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.DG.2.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.DG.2.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.DG.3.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.DG.3.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.dual.0.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.dual.0.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.edge.1.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.edge.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.edge.Full.1.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.edge.Full.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.morley.1.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.morley.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.RT.1.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.RT.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/triangle.RT.2.bas_fun\
	   $(DESTDIR)$(templatedir)/triangle.RT.2.tmp_dof
