define([pkgversion], esyscmd([sh -c "grep Version: DESCRIPTION | cut -d' ' -f2 | tr -d '\n'"]))
AC_INIT([stockassessment],[pkgversion],[an@aqua.dtu.dk])

AC_MSG_NOTICE([${PACKAGE_NAME}: ${PACKAGE_VERSION}])

dnl : ${R_HOME=`R RHOME`}
dnl if test -z "${R_HOME}"; then
dnl   echo "could not determine R_HOME"
dnl   exit 1
dnl fi
dnl CXX=`"${R_HOME}/bin/R" CMD config CXX`
dnl if test -z "$CXX"; then
dnl   AC_MSG_ERROR([No C++ compiler is available])
dnl fi
dnl CXX14=`"${R_HOME}/bin/R" CMD config CXX14`
dnl CXX14PICFLAGS=`"${R_HOME}/bin/R" CMD config CXX14PICFLAGS`
dnl AC_LANG(C++)

dnl AC_MSG_NOTICE([CXX: ${CXX}])
dnl AC_MSG_NOTICE([CXX14: ${CXX14}])
dnl AC_MSG_NOTICE([CXX14PICFLAGS: ${CXX14PICFLAGS}])

dnl if test -z "${CXX14}"; then
dnl    AC_MSG_NOTICE([CXX14 is empty. Trying CXX])
dnl    	save_CXXFLAGS="$CXXFLAGS"
dnl 	CXXFLAGS="$CXXFLAGS -std=c++14 -fpic"
dnl 	AC_MSG_NOTICE([Trying CXXFLAGS: ${CXXFLAGS}])
dnl 	pushdef(prog, [int main(){}])
dnl 	flag_test=0
dnl 	AC_COMPILE_IFELSE([AC_LANG_SOURCE([prog])], [CXX14=$CXX; CXX14PICFLAGS=-fpic],[AC_MSG_FAILURE([Compier with C++14 support is needed])])
dnl 	popdef([prog])
dnl 	CXXFLAGS="$save_CXXFLAGS"
dnl fi

dnl AC_MSG_NOTICE([CXX14: ${CXX14}])
dnl AC_MSG_NOTICE([CXX14PICFLAGS: ${CXX14PICFLAGS}])

AC_ARG_ENABLE(debug,
	AS_HELP_STRING([--enable-debug],[Add debug flags. Note: flags set in a personal Makevars file have precedence]),
	[AS_IF([test $enableval = "no"],[DBG=""],[DBG="-g -O0 -Wall -Wpedantic"]); AC_MSG_NOTICE([Debug flags: ${DBG}])],
	[DBG=""])

AC_ARG_ENABLE(analyze,
	AS_HELP_STRING([--enable-analyze],[Add static analyzer flags. Note: flags set in a personal Makevars file have precedence]),
	[AS_IF([test $enableval = "no"],[ANALYZE=""],[ANALYZE="-fsanitize=address,undefined,leak,shadow-call-stack,pointer-subtract,pointer-compare -fno-sanitize-recover=all -fanalyzer -Wanalyzer-double-fclose -Wanalyzer-double-free -Wanalyzer-exposure-through-output-file -Wanalyzer-file-leak -Wanalyzer-free-of-non-heap -Wanalyzer-malloc-leak -Wanalyzer-mismatching-deallocation -Wanalyzer-null-argument -Wanalyzer-null-dereference -Wanalyzer-possible-null-argument -Wanalyzer-possible-null-dereference -Wanalyzer-shift-count-negative -Wanalyzer-shift-count-overflow -Wanalyzer-stale-setjmp-buffer -Wanalyzer-tainted-array-index -Wanalyzer-too-complex -Wanalyzer-unsafe-call-within-signal-handler -Wanalyzer-use-after-free -Wanalyzer-use-of-pointer-in-stale-stack-frame -Wanalyzer-write-to-const -Wanalyzer-write-to-string-literal"]); AC_MSG_NOTICE([Analyzer flags: ${ANALYZE}])],
	[DBG=""])


AC_ARG_ENABLE(compileunits,
	AS_HELP_STRING([--disable-compileunits],[Compile with header only library instead of using compilation units.]),
	[AS_IF([test $enableval = "yes"],[HEADERONLY="yes"],[HEADERONLY="no"])],
	[HEADERONLY="no"])


AS_IF([test ${HEADERONLY} = "yes"],[AC_MSG_NOTICE([Compiling using header only library])],[AC_MSG_NOTICE([Compiling using compilation units])])


: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
dnl NOTE: There is an issue if TMB is precompiled with CppAD
AC_MSG_NOTICE([${R_HOME}])
PRECOMPILE=`${R_HOME}/bin/Rscript -e "v<-if(file.size(system.file('include/precompile.hpp',package='TMB'))>0){cat('-DTMB_ALREADY_PRECOMPILED')}else{cat('')}"`

AC_MSG_NOTICE([${PRECOMPILE}])

AS_IF([test -z "${PRECOMPILE}"],[AC_MSG_NOTICE([TMB pre-compilation not found])],[AC_MSG_NOTICE([Using TMB pre-compilation])])



AC_SUBST(CXX14)
AC_SUBST(CXX14PICFLAGS)
AC_SUBST(DBG)
AC_SUBST(HEADERONLY)
AC_SUBST(PRECOMPILE)
AC_CONFIG_FILES([src/Makevars])


AC_OUTPUT
