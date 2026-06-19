bin_PROGRAMS = easymesh2mesh\
	easymesh2opendx\
	gmsh2mesh.bin\
	mesh_refine\
	mesh2opendx\
	mesh2easymesh\
	easymesh_refine\
	simpmesh2mesh\
	dbmesh2mesh.bin\
	dxdata2easymeshdata.bin\
	tricube\
	dxdata2tecplotdata.bin\
	hsfc_renum\
	2d_vecfun_dxdata2tecplotdata

easymesh2mesh_SOURCES = easymesh2mesh.cpp
easymesh2opendx_SOURCES = easymesh2opendx.cpp
gmsh2mesh_bin_SOURCES = gmsh2mesh.cpp
mesh_refine_SOURCES = mesh_refine.cpp
mesh2opendx_SOURCES = mesh2opendx.cpp
mesh2easymesh_SOURCES = mesh2easymesh.cpp
easymesh_refine_SOURCES = easymesh_refine.cpp
simpmesh2mesh_SOURCES = simpmesh2mesh.cpp
dbmesh2mesh_bin_SOURCES = dbmesh2mesh.cpp
dxdata2easymeshdata_bin_SOURCES = dxdata2easymeshdata.cpp
tricube_SOURCES = tricube.cpp
dxdata2tecplotdata_bin_SOURCES = dxdata2tecplotdata.cpp
hsfc_renum_SOURCES = hsfc_renum.cpp
2d_vecfun_dxdata2tecplotdata_SOURCES = 2d_vecfun_dxdata2tecplotdata.cpp

bin_CPPFLAGS = -I$(top_srcdir)/include -I$(prefix)/include
LDADD = $(top_builddir)/src/.libs/libAFEPack.so




EXTRA_DIST = gmsh2mesh dbmesh2mesh dxdata2easymeshdata dxdata2tecplotdata

# dist-hook:
# 	cp $(srcdir)/gmsh2mesh $(srcdir)/dbmesh2mesh $(srcdir)/dxdata2easymeshdata $(srcdir)/dxdata2tecplotdata $(distdir)

# install-data-local:
# 	cp $(srcdir)/gmsh2mesh\
# 	   $(srcdir)/dbmesh2mesh\
# 	   $(srcdir)/dxdata2easymeshdata\
# 	   $(srcdir)/dxdata2tecplotdata\
# 	   $(bindir)/
# 	chmod 755 $(bindir)/gmsh2mesh
# 	chmod 755 $(bindir)/dbmesh2mesh
# 	chmod 755 $(bindir)/dxdata2easymeshdata
# 	chmod 755 $(bindir)/dxdata2tecplotdata

# uninstall-local:
# 	rm -f $(bindir)/gmsh2mesh\
#               $(bindir)/dbmesh2mesh\
#               $(bindir)/dxdata2easymeshdata\
# 	      $(bindir)/dxdata2tecplotdata

