
AC_PREREQ(2.61)
AC_INIT(phylosmith, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi

CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

AC_LANG(C++)
AC_REQUIRE_CPP

## Check the C++ compiler using the CXX value set
AC_PROG_CXX
## If it is g++, we have GXX set so let's examine it
if test "${GXX}" = yes; then
    AC_MSG_CHECKING([whether g++ version is sufficient])
    gxx_version=$(${CXX} -v 2>&1 | awk '/^.*g.. version/ {print $3}')
    case ${gxx_version} in
        1.*|2.*|3.*|4.0.*|4.1.*|4.2.*|4.3.*|4.4.*|4.5.*|4.6.*|4.7.0|4.7.1)
             AC_MSG_RESULT([no])
             AC_MSG_WARN([Only g++ version 4.7.2 or greater can be used with RcppArmadillo.])
             AC_MSG_ERROR([Please use a different compiler.])   
        ;;
        4.7.[2-9]|4.8.*|4.9.*|5.*|6.*|7.*|8.*)
             gxx_newer_than_45="-fpermissive"
             AC_MSG_RESULT([(${gxx_version}) yes])
        ;;
    esac
fi

## Default the OpenMP flag to the empty string.
## If and only if OpenMP is found, expand to $(SHLIB_OPENMP_CXXFLAGS)
openmp_flag=""

## Check for broken systems produced by a corporation based in Cupertino
AC_MSG_CHECKING([for macOS])
RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
if test x"${RSysinfoName}" == x"Darwin"; then
   AC_MSG_RESULT([found])
   AC_MSG_WARN([OpenMP unavailable and turned off.])
   have_openmp="#define ARMA_DONT_USE_OPENMP 1"
else
   AC_MSG_RESULT([not found as on ${RSysinfoName}])
   ## Check for OpenMP
   AC_MSG_CHECKING([for OpenMP])
   ## if R has -fopenmp we should be good
   allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
   hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)
   if test x"${hasOpenMP}" == x""; then
      AC_MSG_RESULT([missing])
      have_openmp="#define ARMA_DONT_USE_OPENMP 1"
   else   
      AC_MSG_RESULT([found])
      have_openmp="#define ARMA_USE_OPENMP 1"
      openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
   fi
fi



AC_SUBST([HAVE_OPENMP], ["${have_openmp}"])
AC_SUBST([OPENMP_FLAG], ["${openmp_flag}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT