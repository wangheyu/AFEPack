lib_LTLIBRARIES = libAFEPack_mpi.la
libAFEPack_mpi_la_SOURCES = MPI.cpp\
			MPI_Controller.cpp\
			MPI_LoadBalance.cpp\
			MPI_Migration.cpp

headerincludedir = $(includedir)/AFEPack/mpi
headerinclude_HEADERS = $(top_srcdir)/include/AFEPack/mpi/MPI.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_Controller.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_DOF.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_DOF.templates.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_FaceData.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_HGeometry.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_HGeometry.templates.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_HGeometry_archive.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_LoadBalance.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_MemoryReclaimer.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_Migration.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_PeriodHandler.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_SyncDOF.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_SyncProp.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_UGeometry_archive.h\
		$(top_srcdir)/include/AFEPack/mpi/MPI_ULoadBalance.h


#AM_LDFLAGS = -L/home/gary/local/lib
#AM_CPPFLAGS = -I/home/gary/local/include -I/home/gary/local/include/deal.II -I../include
#AM_CXXFLAGS = -D__SERIALIZATION -DMULTITHREAD -pthread -fno-delete-null-pointer-checks
