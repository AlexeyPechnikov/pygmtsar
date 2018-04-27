dnl -----------------------------------------------------------------
dnl	$Id: configure.ac 228 2015-08-12 00:37:33Z pwessel $
dnl
dnl     This program is free software; you can redistribute it and/or modify
dnl     it under the terms of the GNU General Public License as published by
dnl     the Free Software Foundation; version 2 of the License.
dnl
dnl     This program is distributed in the hope that it will be useful,
dnl     but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl     GNU General Public License for more details.
dnl
dnl     Contact info: David Sandwell, Scripps IO, dsandwell@ucsd.edu
dnl ------------------------------------------------------------------
dnl
dnl Run GNU autoconf to produce a configure script from this input.
dnl ------------------------------------------------------------------
dnl
dnl ------------------------------------------------------------------
dnl Since GMTSAR requires GMT executables version >= 5 we insist GMT
dnl executables must be in the current path; otherwise abort.
dnl 
dnl We also need to gather the following information
dnl 1. Location of GMT include and library files
dnl 2. Location of ORBITS directory
dnl ------------------------------------------------------------------
dnl
AC_INIT(GMTSAR, 5.5.0, gmtsar-l@ucsd.edu)
dnl
dnl ------------------------------------------------------------------
dnl The remainder should need no editing
dnl ------------------------------------------------------------------
dnl
m4_include([ax_lapack.m4])
m4_include([ax_blas.m4])
m4_include([ax_lib_hdf5.m4])
AC_CONFIG_SRCDIR(gmtsar/xcorr.c)
AC_CONFIG_AUX_DIR(`pwd`)
AC_CYGWIN
AC_CANONICAL_SYSTEM
AC_LANG_C
AC_PROG_CC
AC_PROG_CPP
AC_PREFIX_DEFAULT(`pwd`)
AC_PATH_XTRA
dnl
dnl ------------------------------------------------------------------
dnl GMTSAR-specific make macros
dnl ------------------------------------------------------------------
dnl
AC_SUBST(ORBITS_DIR)
AC_SUBST(LALIBS)
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(HDF5_CPPFLAGS)
AC_SUBST(HDF5_LDFLAGS)
AC_SUBST(HDF5_LIBS)
AC_SUBST(LD_OPT)
AC_SUBST(LDFLAGS)
AC_SUBST(GMT_LIB)
AC_SUBST(GMT_INC)
AC_SUBST(TIFF_LIB)
AC_SUBST(TIFF_INC)
AC_SUBST(LIBS)
AC_SUBST(LIBEXT)
AC_SUBST(SL)
AC_SUBST(SL_VERSION)
AC_SUBST(MT)
AC_SUBST(exec_prefix)
AC_SUBST(prefix)
AC_SUBST(rootdir)
AC_SUBST(bindir)
AC_SUBST(GMTSAR_VERSION)
dnl
dnl ------------------------------------------------------------------
dnl Special configure options for GMTSAR installation
dnl ------------------------------------------------------------------
dnl
AC_ARG_WITH(orbits_dir,    [  --with-orbits-dir=      Full path to ORBITS directory[[auto]]], [], [with_orbits_dir=no])
AC_ARG_WITH(tiff_include,  [  --with-tiff-include=    Full path to tiff include directory [[auto]]])
AC_ARG_WITH(tiff_lib,      [  --with-tiff-lib=        Full path to tiff lib directory [[auto]]])
AC_ARG_WITH(gmt_config,    [  --with-gmt-config=      Full path to gmt-config [[auto]]])
AC_ARG_ENABLE(debug,       [  --enable-debug          Compile for debugging instead of optimizing code])
AC_ARG_ENABLE(rpath,	   [  --disable-rpath         Do not hardcode runtime library paths])
AC_ARG_ENABLE(shared,      [  --enable-shared         Build shared (dynamic) GMTSAR libraries instead of static])
AC_ARG_VAR(GMT_INC,Location of GMT headers (compile-time))
AC_ARG_VAR(GMT_LIB,Location of GMT library (compile-time))
dnl
dnl ------------------------------------------------------------------
dnl Check for ar, csh/bash, ld, awk, and ranlib
dnl ------------------------------------------------------------------
dnl
AC_CHECK_PROG(GMT_in_path, gmtvector, 1, 0)
if test "$GMT_in_path" = "X0" ; then
	AC_MSG_RESULT(GMT v5 or later executables are not in your path - aborting)
	exit 1
fi
AC_CHECK_PROGS(PAGER, less more, cat)
AC_CHECK_PROG(AR, ar, ar, :)
AC_CHECK_PROG(CSH, csh, csh, bash, )
AC_CHECK_PROG(LD, ld, ld, :)
AC_CHECK_PROG(GNUTAR, gnutar, gnutar, :)
AC_PROG_LN_S
AC_PROG_INSTALL
AC_PROG_AWK
AC_PROG_RANLIB
dnl
dnl ------------------------------------------------------------------
dnl Get current OS and chip type and set LAPACK
dnl ------------------------------------------------------------------
dnl
os=`uname -s`
AC_MSG_CHECKING(for LAPACK library)
if test "X$os" = "XDarwin" ; then
	MT="-mt"
	osr=`uname -r | cut -d. -f1`
	LALIBS="-framework Accelerate"
else
	MT=""
	AX_LAPACK
	LALIBS="$LAPACK_LIBS $BLAS_LIBS"
fi
AC_MSG_RESULT($LALIBS)
chip=`uname -p`
Sun_rev=`uname -r | grep '^5.'`
rootdir=`pwd`
dnl
dnl ------------------------------------------------------------------
dnl Check for libraries.  Most installs will need to have
dnl -lm but under OS X it is not needed.  Thus, we look
dnl for sin in -lm which should set -lm on most platforms.
dnl ------------------------------------------------------------------
dnl
AC_SEARCH_LIBS(sin, m)
dnl
dnl ------------------------------------------------------------------
dnl Set path to Orbits directory with --with-orbits=dir
dnl ------------------------------------------------------------------
dnl
AC_MSG_CHECKING(for ORBITS directory)
if test "X$with_orbits_dir" = "Xyes"; then
	ORBITS_DIR=$rootdir
	AC_MSG_RESULT([Not specified, default to $ORBITS_DIR])
elif test "X$with_orbits_dir" = "Xno"; then
	ORBITS_DIR=$rootdir
	AC_MSG_RESULT([Not specified, default to $ORBITS_DIR])
else
	ORBITS_DIR=$with_orbits_dir
	AC_MSG_RESULT($ORBITS_DIR)
fi
dnl
dnl ------------------------------------------------------------------
dnl Set the HDF5 include and lib stuff
dnl ------------------------------------------------------------------
AX_LIB_HDF5()
HDF5_LIBS="$HDF5_LIBS -lhdf5_cpp"	# Add the C++ library
AC_MSG_CHECKING(for HDF5 include)
AC_MSG_RESULT($HDF5_CPPFLAGS)
AC_MSG_CHECKING(for HDF5 libs)
AC_MSG_RESULT([$HDF5_LDFLAGS $HDF5_LIBS])
dnl
dnl ------------------------------------------------------------------
dnl Set path to tiff directories with --with-tiff=include|lib
dnl ------------------------------------------------------------------
dnl
AC_MSG_CHECKING(for TIFF include directory)
if test "X$with_TIFF_include" = "X"; then	# Must look for it in standard places
	if [ test -d /sw/include ]; then	# fink
		TIFF_INC=/sw/include
	elif [ test -d /opt/local/include ]; then	# macports
		TIFF_INC=/opt/local/include
	elif [ test -d /usr/local/include ]; then
		TIFF_INC=/usr/local/include
	elif [ test -d /usr/include ]; then	# Ubuntu
		TIFF_INC=/usr/include
	elif [ test -d /opt/csw/include ]; then	# Blastwave
		TIFF_INC=/opt/csw/include
	elif [ test -d /opt/include ]; then
		TIFF_INC=/opt/include
	fi
else
	TIFF_INC=$with_tiff_include
fi
AC_MSG_RESULT($TIFF_INC)
AC_MSG_CHECKING(for TIFF library directory)
if test "X$with_tiff_lib" = "X"; then	# Must look for it in standard places
	if [ test -d /sw/include ]; then	# fink
		TIFF_LIB=/sw/lib
	elif [ test -d /opt/local/include ]; then	# macports
		TIFF_LIB=/opt/local/lib
	elif [ test -d /usr/local/include ]; then
		TIFF_LIB=/usr/local/lib
	elif [ test -d /usr/include ]; then	# Ubuntu
		TIFF_LIB=/usr/lib/x86_64-linux-gnu
	elif [ test -d /opt/csw/include ]; then	# Blastwave
		TIFF_LIB=/opt/csw/lib
	elif [ test -d /opt/include ]; then
		TIFF_LIB=/opt/lib
	fi
else
	TIFF_LIB=$with_tiff_lib/lib
fi
AC_MSG_RESULT($TIFF_LIB)
dnl
dnl ------------------------------------------------------------------
dnl Set path to tiff directories with --with-tiff=include|lib
dnl ------------------------------------------------------------------
dnl
AC_MSG_CHECKING(for TIFF include directory)
if [ test -d /sw/include ]; then	# fink
	TIFF_INC=/sw/include
elif [ test -d /opt/local/include ]; then	# macports
	TIFF_INC=/opt/local/include
elif [ test -d /usr/local/include ]; then
	TIFF_INC=/usr/local/include
elif [ test -d /usr/include ]; then	# Ubuntu
	TIFF_INC=/usr/include
elif [ test -d /opt/csw/include ]; then	# Blastwave
	TIFF_INC=/opt/csw/include
elif [ test -d /opt/include ]; then
	TIFF_INC=/opt/include
fi
AC_MSG_RESULT($TIFF_INC)
AC_MSG_CHECKING(for TIFF library directory)
if test "X$with_tiff_lib" = "X"; then	# Must look for it in standard places
	if [ test -d /sw/include ]; then	# fink
		TIFF_LIB=/sw/lib
	elif [ test -d /opt/local/include ]; then	# macports
		TIFF_LIB=/opt/local/lib
	elif [ test -d /usr/local/include ]; then
		TIFF_LIB=/usr/local/lib
	elif [ test -d /usr/include ]; then	# Ubuntu
		TIFF_LIB=/usr/lib/x86_64-linux-gnu
	elif [ test -d /opt/csw/include ]; then	# Blastwave
		TIFF_LIB=/opt/csw/lib
	elif [ test -d /opt/include ]; then
		TIFF_LIB=/opt/lib
	fi
else
	TIFF_LIB=$with_tiff_lib/lib
fi
AC_MSG_RESULT($TIFF_LIB)
dnl
dnl
dnl ------------------------------------------------------------------
dnl We use gmt-config to set GMT paths and settings
dnl ------------------------------------------------------------------
dnl
if test "X$with_gmt_config" != "Xno" ; then
	AC_PATH_PROG([GMT_CONF], [gmt-config], [$with_gmt_config])
fi
if test -x "$GMT_CONF" && test "X$GMT_INC" = "X" && test "X$GMT_LIB" = "X" ; then
	AC_MSG_CHECKING(for GMT includes)
	GMT_INC=`$GMT_CONF --cflags`
	AC_MSG_RESULT($GMT_INC)
	AC_MSG_CHECKING(for GMT libs)
	GMT_LIB=`$GMT_CONF --libs`
	AC_MSG_RESULT($GMT_LIB)
	GMT_LIB_PATH=`echo $GMT_LIB | sed 's/^-L//;s/\ .*//'`
	rpath="$GMT_LIB_PATH"
	AC_MSG_CHECKING(for 32/64-bit GMT installation)
	bits=`$GMT_CONF --bits`
	GCC_64="-m${bits}"
	LDFLAGS="$LDFLAGS -m${bits}"
	AC_MSG_RESULT(found ${bits}-bit)
fi
dnl
dnl ------------------------------------------------------------------
dnl Allow user to set --enable-debug to pass -g flag
dnl ------------------------------------------------------------------
dnl
AC_MSG_CHECKING(if we should compile with debugging switch -g)
CFLAGS=$ac_save_CFLAGS
CXXFLAGS=$ac_save_CXXFLAGS
if test "X$enable_debug" = "Xyes" ; then
	AC_MSG_RESULT(yes)
	CFLAGS="$CFLAGS -g"
	CXXFLAGS="$CXXFLAGS -g"
else
	CFLAGS=`echo $CFLAGS | sed 's/-g//'`
	CXXFLAGS=`echo $CXXFLAGS | sed 's/-g//'`
	if test "$os" = "Darwin" ; then
		LDFLAGS="$LDFLAGS -Wl,-x"
	else
		LDFLAGS="$LDFLAGS -s"
	fi
	if test "$GCC" = "yes" || test "$os" = "Linux" || test "$os" = "AIX"; then
		CFLAGS="$CFLAGS -O2"
		CXXFLAGS="$CXXFLAGS -O2"
	elif test "$chip" = "alpha" || test "$os" = "AUX" || test "$os" = "HP-UX"; then
		CFLAGS="$CFLAGS"
	else
		CFLAGS="$CFLAGS -O"
		CXXFLAGS="$CXXFLAGS -O"
	fi
	AC_MSG_RESULT(no)
fi
dnl
dnl ------------------------------------------------------------------
dnl Determine LD_OPT for this platform
dnl ------------------------------------------------------------------
dnl
shared_unknown=no
LD_OPT=
LD=ld
SL=so
LIBEXT='$(SL)'
SL_VERSION='$(SL)'
if test "X$enable_shared" = "Xyes" ; then
	AC_MSG_CHECKING(how to make a shared object)
	if test "X$CYGWIN" = "Xyes" ; then
		SL=dll.a
		SL_VERSION=$LIB_MAJOR_VERSION.$SL
		LD='$(CC)'
		LD_OPT="-shared"
	elif test "$os" = "SunOS" ; then
		if test "$Sun_rev" = "" ; then	# SunOS
			LD_OPT="-Bdynamic"
		else				# Solaris
			LD_OPT="-G"
		fi
	elif test "$os" = "Linux" || test "$os" = "GNU/kFreeBSD" ; then	# Need optimization when doing shared
		LD='$(CC)'
		SL_VERSION=$SL.$LIB_MAJOR_VERSION
		LD_OPT='-shared -Wl,-soname=$*.'$SL_VERSION
	elif test "$os" = "FreeBSD" || test "$os" = "NetBSD" ; then
		SL_VERSION=$SL.$LIB_MAJOR_VERSION
		LD_OPT='-shared -Wl,-soname=$*.'$SL_VERSION
	elif test "$os" = "OSF1" ; then
		LD_OPT='-shared -expect_unresolved "*"'
	elif test "$os" = "IRIX64" ; then
		LD_OPT="-shared"
		# Need to remove -L/usr/lib to avoid mixing o32 and n32
		X_LIBS=
	elif test "$os" = "Darwin" ; then	# MacOS X
		SL=dylib
		SL_VERSION=$LIB_MAJOR_VERSION.$SL
		LD='$(CC)'
		LD_OPT='-install_name ${prefix}/lib/$*.'$SL_VERSION' -compatibility_version 4.0.0 -current_version '$LIB_VERSION' -dynamiclib '$OSX_L
		if test $osr -eq 8 ; then	# Tiger needs special option for shared libs
			CFLAGS="$CFLAGS -fno-common"
		fi
	elif test "$os" = "HP-UX" ; then
		LD_OPT="-b"
		SL=sl
	else	# Not yet implemented
		AC_MSG_RESULT(not supported for GMTSAR under $os - building static libs instead)
		enable_shared=no
		shared_unknown=yes
	fi
fi
if test "X$enable_shared" = "Xyes" ; then
	rpath="$rpath"':$(libdir)'
	AC_MSG_RESULT($LD $LD_OPT)
	AC_MSG_CHECKING(the file extension for shared libraries)
	AC_MSG_RESULT($SL)
else
	LIBEXT='a'
fi
dnl
dnl ------------------------------------------------------------------
dnl Determine LDFLAGS for this platform
dnl ------------------------------------------------------------------
dnl
AC_MSG_CHECKING(linker options)
flag=
if test "X$enable_rpath" = "Xno" ; then
	flag=
elif test "$os" = "SunOS" ; then
	if test "$Sun_rev" = "" ; then
		flag="-L"
	else
		flag="-R,"
	fi
elif test "$os" = "Linux" || test "$os" = "FreeBSD" || test "$os" = "NetBSD" || test "$os" = "OSF1" || test "$os" = "IRIX64" ; then
	flag="-rpath,"
elif test "$os" = "HP-UX" ; then
	flag="+b,"
fi
if test "X$flag" != "X" ; then
	LDFLAGS="$LDFLAGS -Wl,$flag$rpath"
fi
AC_MSG_RESULT($LDFLAGS)
dnl
dnl ------------------------------------------------------------------
dnl Determine CFLAGS for this platform
dnl ------------------------------------------------------------------
dnl
AC_MSG_CHECKING(C compiler options for GMTSAR)
if test "$GCC" = "yes"; then		# GNU cc options
	if test "$chip" = "alpha" ; then	# Alpha PC
		CFLAGS="$CFLAGS -mieee -Wall $GCC_64"
		CXXFLAGS="$CXXFLAGS -mieee -Wall $GCC_64"
	else				# MacOSX, Linux, ...
		CFLAGS="$CFLAGS -Wall $GCC_64"
		CXXFLAGS="$CXXFLAGS -Wall $GCC_64"
	fi
	CFLAGS="$CFLAGS -fPIC -fno-strict-aliasing -std=c99"
	CXXFLAGS="$CXXFLAGS -fPIC"
elif test "$os" = "SunOS" ; then	# Sun Workshop CC options
	CFLAGS="$CFLAGS -Xc"
elif test "$os" = "HP-UX" ; then	# HP cc options
	CFLAGS="$CFLAGS -Aa +Ofltacc +O1 +Olibcalls"
	if test "$enable_shared" = "yes" ; then
		CFLAGS="$CFLAGS +z"
	fi
elif test "$os" = "OSF1" ; then		# Dec Alpha OSF cc options
	CFLAGS="$CFLAGS -Olimit 1500 -ieee_with_inexact"
elif test "$os" = "ULTRIX" ; then	# Dec Ultrix cc options
	CFLAGS="$CFLAGS -Olimit 1500"
fi
AC_MSG_RESULT($CFLAGS)
dnl
dnl ------------------------------------------------------------------
dnl Determine CPPFLAGS for this platform
dnl ------------------------------------------------------------------
dnl
AC_MSG_CHECKING(compiler preprocessor flags for GMTSAR)
if test "$os" = "AUX" ; then		# Apple A/UX cc options
	CPPFLAGS="$CPPFLAGS -Daux"
fi
AC_MSG_RESULT($CPPFLAGS)
dnl
dnl ------------------------------------------------------------------
dnl Check if we have GNU make. If not, add CPPFLAGS to CFLAGS
dnl ------------------------------------------------------------------
dnl
AC_MSG_CHECKING(if CPPFLAGS need to be added to CFLAGS)
${MAKE:-make} -p -f /dev/null 2>&1 | grep CPPFLAGS > /dev/null
if test "$?" = "0" ; then
	AC_MSG_RESULT(no)
else
	CFLAGS="$CFLAGS "'$(CPPFLAGS)'
	AC_MSG_RESULT(yes)
fi
GMTSAR_VERSION=AC_PACKAGE_VERSION
dnl
dnl ------------------------------------------------------------------
dnl Update Make macros and preprocess scripts that need the orbits dir
dnl ------------------------------------------------------------------
dnl
AC_OUTPUT(config.mk preproc/ENVI_preproc/scripts/ENVI_pre_process preproc/ENVI_preproc/scripts/ENVI_SLC_pre_process preproc/ERS_preproc/scripts/ERS_pre_process gmtsar/csh/gmtsar_sharedir.csh)
dnl
