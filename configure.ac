#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([AFEPack],[0.0],[rli@math.pku.edu.cn])
AM_INIT_AUTOMAKE([subdir-objects tar-ustar])
LT_INIT
AC_CONFIG_SRCDIR([src/AMGSolver.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])

# Checks for programs.
AC_PROG_CXX([mpicxx])
AC_PROG_CC([mpicc])
AC_PROG_LN_S
dnl  AC_PROG_RANLIB

AC_LANG(C++)

# Checks for libraries.
dnl LIBS="$EXTRA_LIBDIR $LIBS"
dnl Checks for deal.II libraries.
#AC_CHECK_LIB([deal_II], [main],
#  [deal_II_library=yes, LIBS="-ldeal_II $LIBS"],
#  [deal_II_library=no])
# if test "x$deal_II_library" == xno; then
#  AC_CHECK_LIB([base], [main],
#    [deal_II_library=yes, LIBS="-llac -lbase $LIBS"],
#    [deal_II_library=no])
#  if test "x$deal_II_library" == xno; then
#    AC_MSG_ERROR(deal.II library not found)
#  fi
# fi


# Checks for header files.
#AC_CHECK_HEADERS([fcntl.h unistd.h])
#AC_CHECK_HEADERS(deal.II/base/exceptions.h,[],[AC_MSG_ERROR([Could not find deal.II header files])],
#[#define DEBUG
  #include <deal.II/base/exceptions.h>
#])


# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt strchr strrchr strstr])

# preconfig some flags
AC_SUBST([AM_CPPFLAGS],['-I$(top_srcdir)/include'])
AC_SUBST([AM_CXXFLAGS],['-D__SERIALIZATION -DMULTITHREAD -pthread -fopenmp -fno-delete-null-pointer-checks'])


AC_CONFIG_FILES([Makefile
                 tools/Makefile
		 doc/Makefile
		 external/Makefile
		 external/easymesh/Makefile
		 examples/Makefile
		 examples/poisson_equation/Makefile
     examples/poisson_equation_3D/Makefile
     examples/DysonSphere/Makefile
     examples/TaiChiYinYang/Makefile
     examples/elliptic_equation/Makefile
		 examples/step-7/Makefile
		 examples/moving_mesh/Makefile
		 examples/laplacian_evp3D/Makefile
		 examples/coefficient_possion_equation/Makefile
		 examples/local_refine/Makefile
		 examples/local_refine/2d/Makefile
		 examples/local_refine/3d/Makefile
                 include/AFEPack/Makefile
                 include/AFEPack/template/Makefile
                 include/AFEPack/template/four_tetrahedron/Makefile
                 include/AFEPack/template/four_tetrahedron/four_tetrahedron.0/Makefile
                 include/AFEPack/template/four_tetrahedron/four_tetrahedron.1.D/Makefile
                 include/AFEPack/template/four_tetrahedron/four_tetrahedron.1/Makefile
                 include/AFEPack/template/interval/Makefile
                 include/AFEPack/template/interval/interval.1/Makefile
                 include/AFEPack/template/interval/interval.2/Makefile
                 include/AFEPack/template/interval/interval.3/Makefile
                 include/AFEPack/template/interval/interval.4/Makefile
                 include/AFEPack/template/tetrahedron/Makefile
                 include/AFEPack/template/tetrahedron/tetrahedron.0/Makefile
                 include/AFEPack/template/tetrahedron/tetrahedron.1.D/Makefile
                 include/AFEPack/template/tetrahedron/tetrahedron.1/Makefile
                 include/AFEPack/template/tetrahedron/tetrahedron.2.D/Makefile
                 include/AFEPack/template/tetrahedron/tetrahedron.2/Makefile
                 include/AFEPack/template/tetrahedron/tetrahedron.3.D/Makefile
                 include/AFEPack/template/tetrahedron/tetrahedron.edge.1/Makefile
                 include/AFEPack/template/triangle/Makefile
                 include/AFEPack/template/triangle/triangle.0/Makefile
                 include/AFEPack/template/triangle/triangle.1.D/Makefile
                 include/AFEPack/template/triangle/triangle.1/Makefile
                 include/AFEPack/template/triangle/triangle.2/Makefile
                 include/AFEPack/template/triangle/triangle.3/Makefile
                 include/AFEPack/template/triangle/triangle.4/Makefile
                 include/AFEPack/template/triangle/triangle.DG.1/Makefile
                 include/AFEPack/template/triangle/triangle.DG.2/Makefile
                 include/AFEPack/template/triangle/triangle.DG.3/Makefile
                 include/AFEPack/template/triangle/triangle.RT.1/Makefile
                 include/AFEPack/template/triangle/triangle.RT.2/Makefile
                 include/AFEPack/template/triangle/triangle.dual.0/Makefile
                 include/AFEPack/template/triangle/triangle.edge.1/Makefile
                 include/AFEPack/template/triangle/triangle.edge.Full.1/Makefile
                 include/AFEPack/template/triangle/triangle.morley.1/Makefile
                 include/AFEPack/template/twin_tetrahedron/Makefile
                 include/AFEPack/template/twin_tetrahedron/twin_tetrahedron.0/Makefile
                 include/AFEPack/template/twin_tetrahedron/twin_tetrahedron.1.D/Makefile
                 include/AFEPack/template/twin_tetrahedron/twin_tetrahedron.1/Makefile
                 include/AFEPack/template/twin_triangle/Makefile
                 include/AFEPack/template/twin_triangle/twin_triangle.0/Makefile
                 include/AFEPack/template/twin_triangle/twin_triangle.1.D/Makefile
                 include/AFEPack/template/twin_triangle/twin_triangle.1/Makefile
                 include/AFEPack/template/twin_triangle/twin_triangle.2/Makefile
                 include/AFEPack/template/twin_triangle/twin_triangle.DG.1/Makefile
                 include/AFEPack/template/twin_triangle/twin_triangle.RT.1/Makefile
                 include/AFEPack/template/twin_triangle/twin_triangle.morley.1/Makefile
                 include/Makefile
                 src/Makefile
                 src/mpi/Makefile])

AC_OUTPUT
