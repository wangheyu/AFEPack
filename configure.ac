#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([AFEPack],[0.0],[rli@math.pku.edu.cn])
AC_PREFIX_DEFAULT([$HOME/local/AFEPack])
AM_INIT_AUTOMAKE([subdir-objects tar-ustar])
LT_INIT
AC_CONFIG_SRCDIR([src/AMGSolver.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])

# Checks for programs.
AC_PROG_CXX([g++ clang++])
AC_PROG_CC([gcc clang cc])
AC_PROG_LN_S
dnl  AC_PROG_RANLIB

dnl
AC_LANG(C++)

# Check for C++20 support without requiring autoconf-archive macros.
AC_MSG_CHECKING([whether $CXX supports C++20 with -std=c++20])
SAVE_CXX="$CXX"
CXX="$CXX -std=c++20"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <concepts>
template <typename T>
concept Integral = std::integral<T>;]],
    [[static_assert(Integral<int>);
      return 0;]])],
    [AC_MSG_RESULT([yes])
     CXXCPP="$CXXCPP -std=c++20"],
    [AC_MSG_RESULT([no])
     CXX="$SAVE_CXX"
     AC_MSG_WARN([Missing required dependency: a C++ compiler with C++20 support])
     AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install build-essential])
     AC_MSG_ERROR([C++20 support is required])])

# Initialize variables
CBLAS_CPPFLAGS=""
CBLAS_LDFLAGS=""
LAPACKE_CPPFLAGS=""
LAPACKE_LDFLAGS=""
OPENBLAS_LDFLAGS=""
AM_CPPFLAGS=""
AM_LDFLAGS=""
LB_RPATH_FLAGS=""
AFEPACK_APT_CORE="build-essential libopenblas-dev liblapacke-dev libboost-serialization-dev"
AFEPACK_APT_MPI="libopenmpi-dev openmpi-bin"
AFEPACK_APT_GMSH="gmsh"

AC_ARG_ENABLE([mpi],
    [AS_HELP_STRING([--enable-mpi], [build MPI library and install MPI headers @<:@default=no@:>@])],
    [enable_mpi=$enableval],
    [enable_mpi=no])
AS_CASE([$enable_mpi],
    [yes|no], [],
    [AC_MSG_ERROR([--enable-mpi expects yes or no])])

AC_ARG_ENABLE([gmsh-tests],
    [AS_HELP_STRING([--enable-gmsh-tests=auto|yes|no], [run gmsh-based example tests @<:@default=auto@:>@])],
    [enable_gmsh_tests=$enableval],
    [enable_gmsh_tests=auto])
AS_CASE([$enable_gmsh_tests],
    [auto|yes|no], [],
    [AC_MSG_ERROR([--enable-gmsh-tests expects auto, yes, or no])])

# Allow user to specify include path for CBLAS headers
AC_ARG_WITH([cblas-includedir],
    [AS_HELP_STRING([--with-cblas-includedir=DIR], [Specify the directory containing CBLAS header files])],
    [AS_IF([test "x$withval" != xno && test "x$withval" != x],
           [AS_IF([test -d "$withval"],
                  [CBLAS_CPPFLAGS="$withval"],
                  [AC_MSG_WARN([CBLAS include directory $withval does not exist])
                   AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
                   AC_MSG_ERROR([invalid --with-cblas-includedir])])])])

# Allow user to specify library path for general CBLAS libraries
AC_ARG_WITH([cblas-libdir],
    [AS_HELP_STRING([--with-cblas-libdir=DIR], [Specify the directory containing general CBLAS libraries])],
    [AS_IF([test "x$withval" != xno && test "x$withval" != x],
           [AS_IF([test -d "$withval"],
                  [CBLAS_LDFLAGS="$withval"],
                  [AC_MSG_WARN([CBLAS library directory $withval does not exist])
                   AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
                   AC_MSG_ERROR([invalid --with-cblas-libdir])])])])

# Allow user to specify include path for LAPACKE headers
AC_ARG_WITH([lapacke-includedir],
    [AS_HELP_STRING([--with-lapacke-includedir=DIR], [Specify the directory containing LAPACKE header files])],
    [AS_IF([test "x$withval" != xno && test "x$withval" != x],
           [AS_IF([test -d "$withval"],
                  [LAPACKE_CPPFLAGS="$withval"],
                  [AC_MSG_WARN([LAPACKE include directory $withval does not exist])
                   AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
                   AC_MSG_ERROR([invalid --with-lapacke-includedir])])])])

# Allow user to specify library path for general LAPACKE libraries
AC_ARG_WITH([lapacke-libdir],
    [AS_HELP_STRING([--with-lapacke-libdir=DIR], [Specify the directory containing general LAPACKE libraries])],
    [AS_IF([test "x$withval" != xno && test "x$withval" != x],
           [AS_IF([test -d "$withval"],
                  [LAPACKE_LDFLAGS="$withval"],
                  [AC_MSG_WARN([LAPACKE library directory $withval does not exist])
                   AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
                   AC_MSG_ERROR([invalid --with-lapacke-libdir])])])])

# Allow user to specify library path for OpenBLAS libraries
AC_ARG_WITH([openblas-libdir],
    [AS_HELP_STRING([--with-openblas-libdir=DIR], [Specify the directory containing OpenBLAS libraries])],
    [AS_IF([test "x$withval" != xno && test "x$withval" != x],
           [AS_IF([test -d "$withval"],
                  [OPENBLAS_LDFLAGS="$withval"],
                  [AC_MSG_WARN([OpenBLAS library directory $withval does not exist])
                   AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
                   AC_MSG_ERROR([invalid --with-openblas-libdir])])])])

# Apply CBLAS and LAPACKE include paths for header checks
AS_IF([test "x$CBLAS_CPPFLAGS" != x],
    [AM_CPPFLAGS="$AM_CPPFLAGS -I$CBLAS_CPPFLAGS"])
AS_IF([test "x$LAPACKE_CPPFLAGS" != x && test "x$LAPACKE_CPPFLAGS" != "x$CBLAS_CPPFLAGS"],
    [AM_CPPFLAGS="$AM_CPPFLAGS -I$LAPACKE_CPPFLAGS"])
OLD_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $AM_CPPFLAGS"


# Check for CBLAS and LAPACKE headers
AC_CHECK_HEADERS([cblas.h],
    [AC_DEFINE([HAVE_CBLAS_H], [1], [Define if cblas.h is found])],
    [AC_MSG_WARN([Missing required dependency: cblas.h])
     AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
     AC_MSG_WARN([If CBLAS is installed in a custom prefix, pass --with-cblas-includedir=DIR])
     AC_MSG_ERROR([cblas.h not found in paths: $CBLAS_CPPFLAGS])],
    [#include <cblas.h>])

AC_CHECK_HEADERS([lapacke.h],
    [AC_DEFINE([HAVE_LAPACKE_H], [1], [Define if lapacke.h is found])],
    [AC_MSG_WARN([Missing required dependency: lapacke.h])
     AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
     AC_MSG_WARN([If LAPACKE is installed in a custom prefix, pass --with-lapacke-includedir=DIR])
     AC_MSG_ERROR([lapacke.h not found in paths: $LAPACKE_CPPFLAGS])],
    [#include <lapacke.h>])

CPPFLAGS=$OLD_CPPFLAGS

# Apply library paths for library checks
#AM_LDFLAGS="$CBLAS_LDFLAGS $LAPACKE_LDFLAGS $OPENBLAS_LDFLAGS"

# Print the value of MY_VARIABLE
AC_MSG_NOTICE([The value of OPENBLAS_LDFLAGS is: $OPENBLAS_LDFLAGS])
AC_MSG_NOTICE([The value of CBLAS_LDFLAGS is: $CBLAS_LDFLAGS])
AC_MSG_NOTICE([The value of LAPACKE_LDFLAGS is: $LAPACKE_LDFLAGS])
AC_MSG_NOTICE([The value of CBLAS_CPPFLAGS is: $CBLAS_CPPFLAGS])
AC_MSG_NOTICE([The value of LAPACKE_CPPFLAGS is: $LAPACKE_CPPFLAGS])

# Detect the platform
AC_CANONICAL_HOST
case "$host_os" in
  linux*)
    LB_RPATH_FLAGS="-Wl,--enable-new-dtags,-rpath,"
    ;;
  darwin*)
    LB_RPATH_FLAGS="-Wl,-rpath,"
    ;;
  mingw* | cygwin*)
    LB_RPATH_FLAGS="" # Windows typically doesn't use rpath; adjust as needed
    ;;
  freebsd* | openbsd*)
    LB_RPATH_FLAGS="-Wl,-rpath,"
    ;;
  *)
    AC_MSG_ERROR([Unsupported platform: $host_os])
    ;;
esac

# Check for BLAS/LAPACKE.  The OpenBLAS path is a combined-provider mode:
# it must provide both CBLAS and LAPACKE symbols.
SAVE_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $AM_CPPFLAGS"
AS_IF([test "x$OPENBLAS_LDFLAGS" != x],
    [SAVE_LDFLAGS="$LDFLAGS"
     SAVE_LIBS="$LIBS"
     LDFLAGS="$LDFLAGS -L$OPENBLAS_LDFLAGS"
     LIBS="-lopenblas $LIBS"
     AC_MSG_CHECKING([whether OpenBLAS provides CBLAS and LAPACKE symbols])
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>
#include <lapacke.h>]],
        [[double a[1] = {1.0};
          int ipiv[1];
          cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                      1, 1, 1, 1.0, a, 1, a, 1, 0.0, a, 1);
          return LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 1, 1, a, 1, ipiv);]])],
         [
          AC_MSG_RESULT([yes])
       	  AM_LDFLAGS="$AM_LDFLAGS -L$OPENBLAS_LDFLAGS -lopenblas $LB_RPATH_FLAGS$OPENBLAS_LDFLAGS"
          AC_DEFINE([HAVE_OPENBLAS], [1], [Define if OpenBLAS library is found])],
         [AC_MSG_RESULT([no])
          AC_MSG_WARN([Dependency conflict: --with-openblas-libdir was used, but libopenblas in $OPENBLAS_LDFLAGS does not link both cblas_dgemm and LAPACKE_dgetrf])
          AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
          AC_MSG_WARN([Resolution: use an OpenBLAS build with LAPACKE enabled, or omit --with-openblas-libdir and pass matching --with-cblas-libdir/--with-lapacke-libdir paths])
          AC_MSG_ERROR([OpenBLAS/CBLAS/LAPACKE link check failed])])
     LDFLAGS="$SAVE_LDFLAGS"
     LIBS="$SAVE_LIBS"],
    [# Otherwise, check for general CBLAS and LAPACKE libraries.
     # Try libopenblas first (common on Debian/Ubuntu), fall back to libcblas.
     SAVE_LDFLAGS="$LDFLAGS"
     SAVE_LIBS="$LIBS"
     CBLAS_LIBS=""
     AS_IF([test "x$CBLAS_LDFLAGS" != x],
           [LDFLAGS="$LDFLAGS -L$CBLAS_LDFLAGS"])
     LIBS="-lopenblas $LIBS"
     AC_MSG_CHECKING([whether OpenBLAS (system default) provides CBLAS])
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
        [[double a[1] = {1.0};
          cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                      1, 1, 1, 1.0, a, 1, a, 1, 0.0, a, 1);
          return 0;]])],
         [
          AC_MSG_RESULT([yes])
          AS_IF([test "x$CBLAS_LDFLAGS" != x],
                [AM_LDFLAGS="$AM_LDFLAGS -L$CBLAS_LDFLAGS $LB_RPATH_FLAGS$CBLAS_LDFLAGS"])
          CBLAS_LIBS="-lopenblas"
          AC_DEFINE([HAVE_CBLAS], [1], [Define if CBLAS library is found])],
         [AC_MSG_RESULT([no])
          LIBS="-lcblas $SAVE_LIBS"
          AC_MSG_CHECKING([whether general CBLAS library (-lcblas) links])
          AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
             [[double a[1] = {1.0};
               cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                           1, 1, 1, 1.0, a, 1, a, 1, 0.0, a, 1);
               return 0;]])],
              [
               AC_MSG_RESULT([yes])
               AS_IF([test "x$CBLAS_LDFLAGS" != x],
                     [AM_LDFLAGS="$AM_LDFLAGS -L$CBLAS_LDFLAGS $LB_RPATH_FLAGS$CBLAS_LDFLAGS"])
               CBLAS_LIBS="-lcblas"
               AC_DEFINE([HAVE_CBLAS], [1], [Define if CBLAS library is found])],
              [AC_MSG_RESULT([no])
               AC_MSG_WARN([Missing or conflicting dependency: neither -lopenblas nor -lcblas could be linked])
               AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
               AC_MSG_WARN([If CBLAS is installed in a custom prefix, pass --with-cblas-libdir=DIR matching the cblas.h include path])
               AC_MSG_ERROR([CBLAS link check failed])])])
     AM_LDFLAGS="$AM_LDFLAGS $CBLAS_LIBS"
     LDFLAGS="$SAVE_LDFLAGS"
     LIBS="$SAVE_LIBS"

     SAVE_LDFLAGS="$LDFLAGS"
     SAVE_LIBS="$LIBS"
     AS_IF([test "x$LAPACKE_LDFLAGS" != x],
           [LDFLAGS="$LDFLAGS -L$LAPACKE_LDFLAGS"])
     AS_IF([test "x$CBLAS_LDFLAGS" != x],
           [LDFLAGS="$LDFLAGS -L$CBLAS_LDFLAGS"])
     LIBS="-llapacke $CBLAS_LIBS $LIBS"
     AC_MSG_CHECKING([whether LAPACKE library links with $CBLAS_LIBS])
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <lapacke.h>]],
        [[double a[1] = {1.0};
          int ipiv[1];
          return LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 1, 1, a, 1, ipiv);]])],
         [
          AC_MSG_RESULT([yes])
          AS_IF([test "x$LAPACKE_LDFLAGS" != x],
                [AM_LDFLAGS="$AM_LDFLAGS -L$LAPACKE_LDFLAGS $LB_RPATH_FLAGS$LAPACKE_LDFLAGS"])
          AM_LDFLAGS="$AM_LDFLAGS -llapacke"
          AC_DEFINE([HAVE_LAPACKE], [1], [Define if LAPACKE library is found])],
         [AC_MSG_RESULT([no])
          AC_MSG_WARN([Dependency conflict: LAPACKE library could not be linked with $CBLAS_LIBS])
          AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
          AC_MSG_WARN([Resolution: use matching BLAS/LAPACK packages, or pass matching --with-cblas-libdir and --with-lapacke-libdir paths])
          AC_MSG_ERROR([LAPACKE/CBLAS link check failed])])
     LDFLAGS="$SAVE_LDFLAGS"
     LIBS="$SAVE_LIBS"])
CPPFLAGS="$SAVE_CPPFLAGS"

# Check for Boost headers and Boost.Serialization without relying on external
# autoconf-archive Boost macros.
BOOST_CPPFLAGS=""
BOOST_LDFLAGS=""
BOOST_SERIALIZATION_LIB=""
BOOST_PREFIX=""
BOOST_LIBDIR=""

AC_ARG_WITH([boost],
    [AS_HELP_STRING([--with-boost=DIR], [Specify the Boost installation prefix])],
    [AS_IF([test "x$withval" != xno && test "x$withval" != xyes && test "x$withval" != x],
           [AS_IF([test -d "$withval"],
                  [BOOST_PREFIX="$withval"
                   BOOST_CPPFLAGS="-I$withval/include"],
                  [AC_MSG_WARN([Boost prefix $withval does not exist])
                   AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
                   AC_MSG_ERROR([invalid --with-boost])])])])

AC_ARG_WITH([boost-libdir],
    [AS_HELP_STRING([--with-boost-libdir=DIR], [Specify the directory containing Boost libraries])],
    [AS_IF([test "x$withval" != xno && test "x$withval" != x],
           [AS_IF([test -d "$withval"],
                  [BOOST_LIBDIR="$withval"],
                  [AC_MSG_WARN([Boost library directory $withval does not exist])
                   AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
                   AC_MSG_ERROR([invalid --with-boost-libdir])])])])

AS_IF([test "x$BOOST_LIBDIR" = x && test "x$BOOST_PREFIX" != x],
    [for afepack_boost_libdir in "$BOOST_PREFIX/lib" "$BOOST_PREFIX/lib64" "$BOOST_PREFIX/lib/x86_64-linux-gnu"; do
       AS_IF([test -d "$afepack_boost_libdir"],
             [BOOST_LIBDIR="$afepack_boost_libdir"; break])
     done])
AS_IF([test "x$BOOST_LIBDIR" != x],
    [BOOST_LDFLAGS="-L$BOOST_LIBDIR"])

# setup Boost for its LDFLAGS and CPPFLAGS
AM_CPPFLAGS="$AM_CPPFLAGS $BOOST_CPPFLAGS"

SAVE_CPPFLAGS="$CPPFLAGS"
SAVE_LDFLAGS="$LDFLAGS"
SAVE_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
LIBS="-lboost_serialization $LIBS"
AC_MSG_CHECKING([for Boost headers >= 1.65])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <boost/version.hpp>
#if BOOST_VERSION < 106500
#error Boost is too old
#endif]], [[]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_WARN([Missing required dependency: Boost headers >= 1.65])
     AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
     AC_MSG_WARN([If Boost is installed in a custom prefix, pass --with-boost=DIR])
     AC_MSG_ERROR([Boost headers not found or too old])])
AC_CHECK_HEADERS([boost/archive/text_oarchive.hpp],
    [],
    [AC_MSG_WARN([Missing required dependency: Boost.Serialization headers])
     AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
     AC_MSG_ERROR([Boost.Serialization headers not found])])
AC_MSG_CHECKING([whether Boost.Serialization library links])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <boost/archive/text_oarchive.hpp>
#include <sstream>]],
    [[std::ostringstream os;
      boost::archive::text_oarchive oa(os);
      int value = 1;
      oa << value;
      return 0;]])],
    [AC_MSG_RESULT([yes])
     BOOST_SERIALIZATION_LIB="-lboost_serialization"
     AC_SUBST([BOOST_SERIALIZATION_LIB])],
    [AC_MSG_RESULT([no])
     AC_MSG_WARN([Dependency conflict: Boost.Serialization headers were found, but libboost_serialization could not be linked])
     AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_CORE])
     AC_MSG_WARN([Resolution: ensure Boost headers and libboost_serialization are from the same installation prefix; use --with-boost=DIR and --with-boost-libdir=DIR if needed])
     AC_MSG_ERROR([Boost.Serialization link check failed])])
CPPFLAGS="$SAVE_CPPFLAGS"
LDFLAGS="$SAVE_LDFLAGS"
LIBS="$SAVE_LIBS"
BOOST_RPATH=$(echo "$BOOST_LDFLAGS" | sed 's/-L\([^ ]*\)/-Wl,-rpath,\1/g')
AM_LDFLAGS="$AM_LDFLAGS $BOOST_LDFLAGS $BOOST_SERIALIZATION_LIB $BOOST_RPATH"
AC_SUBST([BOOST_CPPFLAGS])
AC_SUBST([BOOST_LDFLAGS])

AC_MSG_NOTICE([The value of BOOST_LDFLAGS is: $BOOST_LDFLAGS])
AC_MSG_NOTICE([The value of BOOST_CPPFLAGS is: $BOOST_CPPFLAGS])

AC_CHECK_PROG([GMSH], [gmsh], [gmsh], [])
AS_IF([test "x$enable_gmsh_tests" = xyes && test "x$GMSH" = x],
    [AC_MSG_WARN([Missing test dependency: gmsh])
     AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_GMSH])
     AC_MSG_ERROR([--enable-gmsh-tests=yes requires gmsh])])
AS_IF([test "x$enable_gmsh_tests" = xno],
    [GMSH_TESTS=no],
    [AS_IF([test "x$GMSH" = x],
        [GMSH_TESTS=no
         AC_MSG_WARN([gmsh not found; 3D example run tests will be skipped])
         AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_GMSH])],
        [GMSH_TESTS=yes])])
AC_SUBST([GMSH])
AC_SUBST([GMSH_TESTS])

MPICXX=""
AS_IF([test "x$enable_mpi" = xyes],
    [AC_CHECK_PROGS([MPICXX], [mpicxx mpic++ mpiCC], [])
     AS_IF([test "x$MPICXX" = x],
        [AC_MSG_WARN([Missing MPI dependency: mpicxx])
         AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_MPI])
         AC_MSG_ERROR([--enable-mpi requires an MPI C++ compiler])])
     AC_MSG_CHECKING([whether $MPICXX can compile MPI code with C++20])
     cat > conftest_mpi.cpp <<'_ACEOF'
#include <mpi.h>
#include <concepts>
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Finalize();
  return 0;
}
_ACEOF
     AS_IF([$MPICXX -std=c++20 -c conftest_mpi.cpp >/dev/null 2>&1],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_WARN([Dependency conflict: $MPICXX cannot compile a C++20 MPI program])
         AC_MSG_WARN([Debian/Ubuntu install suggestion: sudo apt install $AFEPACK_APT_MPI])
         AC_MSG_WARN([Resolution: install an MPI implementation with C++ compiler wrappers, or set MPICXX=/path/to/mpicxx])
         rm -f conftest_mpi.cpp conftest_mpi.o
         AC_MSG_ERROR([MPI C++20 compile check failed])])
     rm -f conftest_mpi.cpp conftest_mpi.o])
AC_SUBST([MPICXX])
AM_CONDITIONAL([ENABLE_MPI], [test "x$enable_mpi" = xyes])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt strchr strrchr strstr])


# Detect the platform
AC_CANONICAL_HOST
case "$host_os" in
  linux*)
    RPATH_FLAGS="-Wl,--enable-new-dtags,-rpath,${libdir}"
    ;;
  darwin*)
    RPATH_FLAGS="-Wl,-rpath,@loader_path/../lib"
    ;;
  mingw* | cygwin*)
    RPATH_FLAGS="" # Windows typically doesn't use rpath; adjust as needed
    ;;
  freebsd* | openbsd*)
    RPATH_FLAGS="-Wl,-rpath,\$$ORIGIN/../lib"
    ;;
  *)
    AC_MSG_ERROR([Unsupported platform: $host_os])
    ;;
esac
AC_SUBST(RPATH_FLAGS)

# Add --enable-debug option
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable debugging, disables optimizations])],
    [enable_debug=$enableval],
    [enable_debug=no])

# Set compiler and linker flags
# Clear CXXFLAGS to avoid default -O2 -g
CXXFLAGS=""
AS_IF([test "x$enable_debug" = "xyes"],
    [
        AM_CXXFLAGS="-O0 -g -D__SERIALIZATION -DMULTITHREAD -pthread -std=c++20 -Wall -Wextra"
        AC_MSG_NOTICE([Configuring for debug mode with -g -O0])
    ],
    [
        AM_CXXFLAGS="-O2 -D__SERIALIZATION -DMULTITHREAD -pthread -std=c++20 -Wall"
        AC_MSG_NOTICE([Configuring for release mode with -O2])
    ])


# Export variables to Makefile
AC_SUBST([AM_LDFLAGS])
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(top_srcdir)/include"
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_CXXFLAGS])

AC_CONFIG_FILES([examples/check-2d-examples.sh], [chmod +x examples/check-2d-examples.sh])
AC_CONFIG_FILES([examples/check-3d-examples.sh], [chmod +x examples/check-3d-examples.sh])

AC_CONFIG_FILES([Makefile
                 tools/Makefile
		 doc/Makefile
		 external/Makefile
		 external/easymesh/Makefile
		 examples/Makefile
		 examples/poisson_equation/Makefile
		 examples/poisson_equation/Makefile.sample
		 examples/poisson_equation/run.sh.sample		 		 
		 examples/poisson_equation_3D/Makefile
		 examples/poisson_equation_3D/Makefile.sample
		 examples/poisson_equation_3D/run.sh.sample		 		 
		 examples/elliptic_equation/Makefile
		 examples/elliptic_equation/Makefile.sample
		 examples/elliptic_equation/run.sh.sample
		 examples/TaiChiYinYang/Makefile
		 examples/TaiChiYinYang/Makefile.sample
		 examples/TaiChiYinYang/run.sh.sample		 
		 examples/DysonSphere/Makefile
		 examples/DysonSphere/Makefile.sample
		 examples/DysonSphere/run.sh.sample	 
		 tools/gmsh2mesh.sh
                 include/AFEPack/Makefile
		 include/AFEPack/mpi/Makefile
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
