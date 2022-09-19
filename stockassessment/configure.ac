define([pkgversion], esyscmd([sh -c "grep Version: DESCRIPTION | cut -d' ' -f2 | tr -d '\n'"]))
AC_INIT([stockassessment],[pkgversion],[an@aqua.dtu.dk])

AC_MSG_NOTICE([${PACKAGE_NAME}: ${PACKAGE_VERSION}])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
if test -z "$CXX"; then
  AC_MSG_ERROR([No C++ compiler is available])
fi
CXX14=`"${R_HOME}/bin/R" CMD config CXX14`
CXX14PICFLAGS=`"${R_HOME}/bin/R" CMD config CXX14PICFLAGS`
AC_LANG(C++)

AC_MSG_NOTICE([CXX: ${CXX}])
AC_MSG_NOTICE([CXX14: ${CXX14}])
AC_MSG_NOTICE([CXX14PICFLAGS: ${CXX14PICFLAGS}])

if test -z "${CXX14}"; then
   AC_MSG_NOTICE([CXX14 is empty. Trying CXX])
   	save_CXXFLAGS="$CXXFLAGS"
	CXXFLAGS="$CXXFLAGS -std=c++14 -fpic"
	AC_MSG_NOTICE([Trying CXXFLAGS: ${CXXFLAGS}])
	pushdef(prog, [int main(){}])
	flag_test=0
	AC_COMPILE_IFELSE([AC_LANG_SOURCE([prog])], [CXX14=$CXX; CXX14PICFLAGS=-fpic],[AC_MSG_FAILURE([Compier with C++14 support is needed])])
	popdef([prog])
	CXXFLAGS="$save_CXXFLAGS"
fi

AC_MSG_NOTICE([CXX14: ${CXX14}])
AC_MSG_NOTICE([CXX14PICFLAGS: ${CXX14PICFLAGS}])


AC_SUBST(CXX14)
AC_SUBST(CXX14PICFLAGS)
AC_CONFIG_FILES([src/Makevars])


AC_OUTPUT
