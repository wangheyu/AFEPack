bin_PROGRAMS =  easymesh2mesh\
	        easymesh2opendx\
	        gmsh2mesh\
		mesh_refine\
		mesh2opendx\
		mesh2easymesh\
		easymesh_refine\
		simpmesh2mesh\
		dbmesh2mesh\
		dxdata2easymeshdata\
		tricube\
		dxdata2tecplotdata\
		hsfc_renum\
		2d_vecfun_dxdata2tecplotdata

EXTRA_DIST = gmsh2mesh.sh.in

nodist_bin_SCRIPTS = gmsh2mesh.sh

easymesh2mesh_SOURCES = easymesh2mesh.cpp
easymesh2opendx_SOURCES = easymesh2opendx.cpp
gmsh2mesh_SOURCES = gmsh2mesh.cpp
mesh_refine_SOURCES = mesh_refine.cpp
mesh2opendx_SOURCES = mesh2opendx.cpp
mesh2easymesh_SOURCES = mesh2easymesh.cpp
easymesh_refine_SOURCES = easymesh_refine.cpp
simpmesh2mesh_SOURCES = simpmesh2mesh.cpp
dbmesh2mesh_SOURCES = dbmesh2mesh.cpp
dxdata2easymeshdata_SOURCES = dxdata2easymeshdata.cpp
tricube_SOURCES = tricube.cpp
dxdata2tecplotdata_SOURCES = dxdata2tecplotdata.cpp
hsfc_renum_SOURCES = hsfc_renum.cpp
2d_vecfun_dxdata2tecplotdata_SOURCES = 2d_vecfun_dxdata2tecplotdata.cpp

LDADD =$(top_builddir)/src/libAFEPack.la
AM_LDFLAGS+=$(RPATH_FLAGS)
