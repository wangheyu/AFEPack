templatedir = $(includedir)/AFEPack/template/interval

template_LTLIBRARIES = libinterval.geometry.la

libinterval_geometry_la_SOURCES = interval.crd_trs.c\
	interval.to2d.crd_trs.c\
	interval.volume.c

SUBDIRS = interval.1\
	interval.2\
	interval.3\
	interval.4


EXTRA_DIST = interval.crd_trs\
	interval.tmp_geo\
	interval.to2d.crd_trs

install-data-local:
	cp $(srcdir)/interval.crd_trs\
	   $(srcdir)/interval.tmp_geo\
	   $(srcdir)/interval.to2d.crd_trs\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/interval.1/interval.1.bas_fun\
	   $(srcdir)/interval.1/interval.1.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/interval.2/interval.2.bas_fun\
	   $(srcdir)/interval.2/interval.2.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/interval.3/interval.3.bas_fun\
	   $(srcdir)/interval.3/interval.3.tmp_dof\
	   $(DESTDIR)$(templatedir)/
	cp $(srcdir)/interval.4/interval.4.bas_fun\
	   $(srcdir)/interval.4/interval.4.tmp_dof\
	   $(DESTDIR)$(templatedir)/

uninstall-local:
	rm -f $(DESTDIR)$(templatedir)/interval.crd_trs\
              $(DESTDIR)$(templatedir)/interval.to2d.crd_trs\
              $(DESTDIR)$(templatedir)/interval.tmp_geo
	rm -f $(DESTDIR)$(templatedir)/interval.1.bas_fun\
	      $(DESTDIR)$(templatedir)/interval.1.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/interval.2.bas_fun\
	      $(DESTDIR)$(templatedir)/interval.2.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/interval.3.bas_fun\
	      $(DESTDIR)$(templatedir)/interval.3.tmp_dof
	rm -f $(DESTDIR)$(templatedir)/interval.4.bas_fun\
	      $(DESTDIR)$(templatedir)/interval.4.tmp_dof




#filedir = $(top_srcdir)/library/include/AFEPack/template/interval



# dist-hook:
# 	cp $(srcdir)/interval.crd_trs $(distdir)
# 	cp $(srcdir)/interval.tmp_geo $(distdir)
# 	cp $(srcdir)/interval.to2d.crd_trs $(distdir)


# install-data-local:
# 	cp $(filedir)/interval.crd_trs\
# 	   $(filedir)/interval.tmp_geo\
# 	   $(filedir)/interval.to2d.crd_trs\
# 	   $(templatedir)/
# 	cp $(filedir)/interval.1/interval.1.bas_fun\
# 	   $(filedir)/interval.1/interval.1.tmp_dof\
# 	   $(templatedir)/
# 	cp $(filedir)/interval.2/interval.2.bas_fun\
# 	   $(filedir)/interval.2/interval.2.tmp_dof\
# 	   $(templatedir)/
# 	cp $(filedir)/interval.3/interval.3.bas_fun\
# 	   $(filedir)/interval.3/interval.3.tmp_dof\
# 	   $(templatedir)/
# 	cp $(filedir)/interval.4/interval.4.bas_fun\
# 	   $(filedir)/interval.4/interval.4.tmp_dof\
# 	   $(templatedir)/

# uninstall-local:
# 	rm -f $(templatedir)/interval.crd_trs\
#               $(templatedir)/interval.to2d.crd_trs\
#               $(templatedir)/interval.tmp_geo
# 	rm -f $(templatedir)/interval.1.bas_fun\
# 	      $(templatedir)/interval.1.tmp_dof
# 	rm -f $(templatedir)/interval.2.bas_fun\
# 	      $(templatedir)/interval.2.tmp_dof
# 	rm -f $(templatedir)/interval.3.bas_fun\
# 	      $(templatedir)/interval.3.tmp_dof
# 	rm -f $(templatedir)/interval.4.bas_fun\
# 	      $(templatedir)/interval.4.tmp_dof

