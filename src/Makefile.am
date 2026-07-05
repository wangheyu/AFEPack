lib_LTLIBRARIES = libAFEPack.la
libAFEPack_la_SOURCES = $(srcdir)/AMGSolver.cpp\
			$(srcdir)/BilinearOperator.cpp\
			$(srcdir)/BoundaryCondition.cpp\
			$(srcdir)/DBMesh.cpp\
			$(srcdir)/DGFEMSpace.cpp\
			$(srcdir)/EasyMesh.cpp\
			$(srcdir)/FEMSpace.cpp\
			$(srcdir)/FEMSpace.vector_value.1d.cpp\
			$(srcdir)/FEMSpace.vector_value.2d.cpp\
			$(srcdir)/FEMSpace.vector_value.3d.cpp\
			$(srcdir)/Functional.cpp\
			$(srcdir)/Geometry.1d.cpp\
			$(srcdir)/Geometry.2d.cpp\
			$(srcdir)/Geometry.3d.cpp\
			$(srcdir)/Geometry.cpp\
			$(srcdir)/GmshMesh.cpp\
			$(srcdir)/HGeometry.1d.cpp\
			$(srcdir)/HGeometry.2d.cpp\
			$(srcdir)/HGeometry.3d.cpp\
			$(srcdir)/Migration.cpp\
			$(srcdir)/Miscellaneous.cpp\
			$(srcdir)/MovingMesh2D.cpp\
			$(srcdir)/MovingMesh3D.cpp\
			$(srcdir)/MovingMeshFB.cpp\
			$(srcdir)/MovingMeshFB3D.cpp\
			$(srcdir)/Operator.discretize.cpp\
			$(srcdir)/Operator.interpolate.cpp\
			$(srcdir)/Operator.project.cpp\
			$(srcdir)/PropertyTable.cpp\
			$(srcdir)/SparseMatrixTool.cpp\
			$(srcdir)/TemplateElement.cpp\
			$(srcdir)/TemplateElement.vector_value.l1.cpp\
			$(srcdir)/TemplateElement.vector_value.l2.cpp\
			$(srcdir)/TemplateElement.vector_value.l3.cpp\
			$(srcdir)/Thread.cpp

# libAFEPack_la_CPPFLAGS = -I$(top_srcdir)/include

EXTRA_DIST = $(srcdir)/HGeometry.cpptemplate

if ENABLE_MPI
SUBDIRS = mpi
else
SUBDIRS =
endif

DIST_SUBDIRS = mpi
