lib_LTLIBRARIES = libAFEPack.la
libAFEPack_la_LIBADD = -lopenblas
libAFEPack_la_SOURCES =AMGSolver.cpp\
			BilinearOperator.cpp\
			BoundaryCondition.cpp\
			DBMesh.cpp\
			DGFEMSpace.cpp\
			DenseMatrix.cpp\
			EasyMesh.cpp\
			FEMSpace.cpp\
			FEMSpace.vector_value.1d.cpp\
			FEMSpace.vector_value.2d.cpp\
			FEMSpace.vector_value.3d.cpp\
			Functional.cpp\
			Geometry.1d.cpp\
			Geometry.2d.cpp\
			Geometry.3d.cpp\
			Geometry.cpp\
			GmshMesh.cpp\
			HGeometry.1d.cpp\
			HGeometry.2d.cpp\
			HGeometry.3d.cpp\
			LOBPCG.cpp\
			Migration.cpp\
			Miscellaneous.cpp\
			MovingMesh2D.cpp\
			MovingMesh3D.cpp\
			MovingMeshFB.cpp\
			MovingMeshFB3D.cpp\
			Operator.discretize.cpp\
			Operator.interpolate.cpp\
			Operator.project.cpp\
			PropertyTable.cpp\
			SparseMatrix.cpp\
			SparseMatrixTool.cpp\
			SparsityPattern.cpp\
			TemplateElement.cpp\
			TemplateElement.vector_value.l1.cpp\
			TemplateElement.vector_value.l2.cpp\
			TemplateElement.vector_value.l3.cpp\
			Thread.cpp\
			Vector.cpp

EXTRA_DIST = HGeometry.cpptemplate


# install-exec-hook:
# 	$(LN_S) .libs/libAFEPack.so $(top_srcdir)/lib/

SUBDIRS = mpi
