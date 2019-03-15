dnl Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
AC_INIT(AFEPack, 1.8, [rli@math.pku.edu.cn])
AC_CONFIG_SRCDIR([library])
AM_CONFIG_HEADER(config.h)
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE

AFEPACK_PATH=`pwd`
AC_SUBST(AFEPACK_PATH)

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_LN_S
AC_CHECK_TOOL(AR, ar, :)
AC_CHECK_TOOL(DOCGEN, doxygen, :)
AC_CHECK_TOOL(PERL, perl, :)
AC_CHECK_TOOL(SED, sed, :)
AC_CHECK_PROGS(MAKE, [gmake make])
AC_LANG(C++)

dnl libdl
AC_TRY_CPP([#include <dlfcn.h>],
  dl_header=yes, dl_header=no)
AC_MSG_CHECKING([for dlfcn.h])
if test "x$dl_header" != xno; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_ERROR([dlfcn.h not found])
fi

AC_CHECK_LIB([dl], [dlopen], [], [AC_ERROR([library -ldl not found])])
dnl ----------------------

dnl lib math
AC_CHECK_LIB([m], [sin], [], [AC_ERROR([library -lm not found])])
dnl ----------------------

THIS_FLAG="$EXTRA_INCDIR"
CFLAGS="$THIS_FLAG $CFLAGS"
CXXFLAGS="$THIS_FLAG $CXXFLAGS"
CPPFLAGS="$THIS_FLAG $CPPFLAGS"

enable_shared=yes
AC_SUBST(enable_shared)
LIBS="$EXTRA_LIBS $LIBS"

dnl multi-thread option
AC_ARG_ENABLE(thread,
  AC_HELP_STRING([--enable-thread], [enable POSIX thread support(default: yes)]),,
  [enable_thread=yes])
AC_MSG_CHECKING([if thread support was requested])
if test "x$enable_thread" != xno; then
  AC_MSG_RESULT(yes)
  THIS_FLAG="-D__SERIALIZATION__ -DMULTITHREAD -pthread -fPIC"
  CFLAGS="$CFLAGS $THIS_FLAG"
  CXXFLAGS="$CXXFLAGS $THIS_FLAG -std=c++11"
  CPPFLAGS="$CPPFLAGS $THIS_FLAG"
else
  AC_MSG_RESULT(no)
fi
AC_SUBST(enable_thread)

dnl MPI library
AC_MSG_CHECKING([for MPI])
AC_ARG_WITH([mpi], AC_HELP_STRING([--with-mpi={yes|no|PATH}], [Specify if MPI is enabled or path to mpi installation]), [], [with_mpi=check])
if test "x$with_mpi" = xno; then
  with_mpi=check
fi
if test "x$with_mpi" != xcheck; then
  if test "$with_mpi" = yes; then
    MPICC=`which mpicc`
    MPICXX=`which mpicxx`
  else
    MPICC="$with_mpi/bin/mpicc"
    MPICXX="$with_mpi/bin/mpicxx"
  fi
  CXX_backup=$CXX
  CXX=$MPICXX
  AC_TRY_LINK([#include <mpi.h>],
    [MPI_Init(NULL, NULL); MPI_Finalize()], mpi_lib=yes, mpi_lib=no)
  CXX=$CXX_backup
  if test "x$mpi_lib" != xno; then
    AC_MSG_RESULT(yes)
  else
   AC_MSG_RESULT(no)
   AC_ERROR([MPI not found])
  fi
else
  AC_MSG_RESULT(no)
  MPICC=""
  MPICXX=""
fi
AC_SUBST(MPICC)
AC_SUBST(MPICXX)
AC_SUBST(with_mpi)
dnl ----------------------

dnl boost headers
AC_ARG_WITH([boost], AC_HELP_STRING([--with-boost=PATH], [Specify path to boost installation]), [], [with_boost=check])
if test "x$with_boost" != xcheck; then
  THIS_FLAG="-I$with_boost/include"
  CFLAGS="$THIS_FLAG $CFLAGS"
  CXXFLAGS="$THIS_FLAG $CXXFLAGS"
  CPPFLAGS="$THIS_FLAG $CPPFLAGS"
fi
AC_TRY_CPP([#include <boost/serialization/array.hpp>],
  boost_header=yes, boost_header=no)
AC_MSG_CHECKING([for boost headers])
if test "x$dl_header" != xno; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_ERROR([boost headers not found])
fi
AC_SUBST(with_boost)
dnl ----------------------

dnl Header for deal.II
AC_ARG_WITH([dealii], AC_HELP_STRING([--with-dealii=PATH], [Specify path to deal.II installation]), [], [with_dealii=check])
if test "x$with_dealii" != xcheck; then
  THIS_FLAG="-I$with_dealii/base/include -I$with_dealii/lac/include -I$with_dealii/include -I$with_dealii/include/deal.II"
  CFLAGS="$THIS_FLAG $CFLAGS"
  CXXFLAGS="$THIS_FLAG $CXXFLAGS"
  CPPFLAGS="$THIS_FLAG $CPPFLAGS"
fi
AC_TRY_CPP([#include <base/config.h>],
  deal_II_header=yes, deal_II_header=no)
AC_MSG_CHECKING([for deal.II header files])
if test "x$deal_II_header" != xno; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_ERROR([deal.II header not found])
fi

dnl Libraries for deal.II
AC_MSG_CHECKING([for deal.II library (-ldeal_II)])
LIBS="-L$with_dealii/lib $LIBS"
LIBS_backup=$LIBS
LIBS="-ldeal_II $LIBS_backup"
AC_TRY_LINK([#include <base/config.h>], [],
  deal_II_library=yes, deal_II_library=no)
if test "x$deal_II_library" != xno; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(not found)
  AC_MSG_CHECKING([for deal.II library (-llac -lbase -ltbb)])
  LIBS="-llac -lbase -ltbb $LIBS_backup"
  AC_TRY_LINK([#include <base/config.h>], [],
  deal_II_library=yes, deal_II_library=no)
  if test "x$deal_II_library" != xno; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(not found)
    AC_ERROR([deal.II library not found])
  fi
fi
dnl -----------------------

THIS_FLAG="-ftemplate-depth-256 -I$AFEPACK_PATH/library/include"
CFLAGS="$CFLAGS $THIS_FLAG"
CXXFLAGS="$CXXFLAGS $THIS_FLAG"
CPPFLAGS="$CPPFLAGS $THIS_FLAG"
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(CPPFLAGS)

AC_SUBST(LIBS)

dnl Checks for header files.
AC_HEADER_STDC

AC_OUTPUT(Makefile \
  template/Makefile \
  example/Makefile \
  example/local_refine/Makefile \
  contrib/Makefile \
  Make.global_options)

