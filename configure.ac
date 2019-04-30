
AC_PREREQ(2.69)
AC_INIT(phylosmith)
AC_CONFIG_SRCDIR([src])
AC_LANG(C++)
AC_REQUIRE_CPP

: ${R_HOME=$(R RHOME)}
AC_MSG_CHECKING([R_HOME already set])
if test -z "${R_HOME}"; then
    AC_MSG_RESULT([No])
    AC_MSG_ERROR([Cannot proceed.])
fi
AC_MSG_RESULT([${R_HOME}])
R_BIN="${R_HOME}/bin/R"

AC_MSG_CHECKING([R CXX11])
if test -z "$CXX11"; then
  CXX11=`"${R_HOME}/bin/R" CMD config CXX11`
  AC_MSG_RESULT([Setting ${CXX11}])
else
  AC_MSG_RESULT([Using pre-set ${CXX11}])
fi

AC_MSG_CHECKING([R CXX11STD])
if test -z "$CXX11STD"; then
  CXX11STD=`"${R_HOME}/bin/R" CMD config CXX11STD`
  AC_MSG_RESULT([Setting R ${CXX11STD}])
else
  AC_MSG_RESULT([Using existing ${CXX11STD}])
fi

AC_MSG_CHECKING([R CXX11FLAGS])
if test -z "$CXX11FLAGS"; then
  CXX11FLAGS=`"${R_HOME}/bin/R" CMD config CXX11FLAGS`
  AC_MSG_RESULT([Setting R ${CXX11FLAGS}])
else
  AC_MSG_RESULT([Using existing ${CXX11FLAGS}])
fi

CXX1X=`"${R_HOME}/bin/R" CMD config CXX1X`
CXX1XSTD=`"${R_HOME}/bin/R" CMD config CXX1XSTD`
CXX="${CXX1X} ${CXX1XSTD}"
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX1XFLAGS`

LIBS="${PKG_LIBS}"

OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
AC_ARG_ENABLE([openmp],
[AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
if test "$enable_openmp" != no; then
AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
  [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
  [AC_LINK_IFELSE([_AC_LANG_OPENMP],
 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
  for ac_option in -fopenmp -xopenmp -qopenmp \
                       -openmp -mp -omp -qsmp=omp -homp \
		   -fopenmp=libomp \
                       -Popenmp --openmp; do
    ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
    _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
    AC_LINK_IFELSE([_AC_LANG_OPENMP],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
    _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
    if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
      break
    fi
  done])])
case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
  "none needed" | unsupported)
;; #(
  *)
OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp ;;
esac
fi

AC_SUBST(OPENMP_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT