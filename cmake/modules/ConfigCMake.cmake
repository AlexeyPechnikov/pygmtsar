#
# $Id$
#
# Useful CMake variables.
#
# There are three configuration files:
#   1) "ConfigDefault.cmake" - is version controlled and used to add new default
#      variables and set defaults for everyone.
#   2) "ConfigUser.cmake" in the source tree - is not version controlled
#      (currently listed in svn:ignore property) and used to override defaults on
#      a per-user basis.
#   3) "ConfigUser.cmake" in the build tree - is used to override
#      "ConfigUser.cmake" in the source tree.
#
# NOTE: If you want to change CMake behaviour just for yourself then copy
#      "ConfigUserTemplate.cmake" to "ConfigUser.cmake" and then edit
#      "ConfigUser.cmake" (not "ConfigDefault.cmake" or "ConfigUserTemplate.cmake").
#
include ("${CMAKE_SOURCE_DIR}/cmake/ConfigDefault.cmake")

# If "ConfigUser.cmake" doesn't exist then create one for convenience.
if (EXISTS "${CMAKE_SOURCE_DIR}/cmake/ConfigUser.cmake")
	include ("${CMAKE_SOURCE_DIR}/cmake/ConfigUser.cmake")
endif (EXISTS "${CMAKE_SOURCE_DIR}/cmake/ConfigUser.cmake")

# If you've got a 'ConfigUser.cmake' in the build tree then that overrides the
# one in the source tree.
if (EXISTS "${CMAKE_BINARY_DIR}/cmake/ConfigUser.cmake")
	include ("${CMAKE_BINARY_DIR}/cmake/ConfigUser.cmake")
endif (EXISTS "${CMAKE_BINARY_DIR}/cmake/ConfigUser.cmake")

###########################################################
# Do any needed processing of the configuration variables #
###########################################################

# Build type
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

# Here we change it to add the SVN revision number for non-public releases - see Package.cmake for
# why this has to be done here.
set (GMTSAR_PACKAGE_VERSION_WITH_SVN_REVISION ${GMTSAR_PACKAGE_VERSION})
# Add the Subversion version number to the package filename if this is a non-public release.
# A non-public release has an empty 'GMTSAR_SOURCE_CODE_CONTROL_VERSION_STRING' variable in 'ConfigDefault.cmake'.
set (HAVE_SVN_VERSION)
if (NOT GMTSAR_SOURCE_CODE_CONTROL_VERSION_STRING)
	# Get the location, inside the staging area location, to copy the application bundle to.
	execute_process (
		COMMAND svnversion ${GMTSAR_SOURCE_DIR}
		RESULT_VARIABLE SVN_VERSION_RESULT
		OUTPUT_VARIABLE SVN_VERSION_OUTPUT
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	if (SVN_VERSION_RESULT)
		message (STATUS "Unable to determine svn version number for non-public release - ignoring.")
	else (SVN_VERSION_RESULT)
		if (SVN_VERSION_OUTPUT MATCHES "Unversioned")
			message (STATUS "Unversioned source tree, non-public release.")
		else (SVN_VERSION_OUTPUT MATCHES "Unversioned")
			# The 'svnversion' command can output a range of revisions with a colon
			# separator - but this causes problems with filenames so we'll remove the
			# colon and the end revision after it.
			string (REGEX REPLACE ":.*$" "" SVN_VERSION ${SVN_VERSION_OUTPUT})
			if (NOT SVN_VERSION STREQUAL exported)
				# Set the updated package version.
				set (GMTSAR_PACKAGE_VERSION_WITH_SVN_REVISION "${GMTSAR_PACKAGE_VERSION}_r${SVN_VERSION}")
				set (HAVE_SVN_VERSION TRUE)
			endif (NOT SVN_VERSION STREQUAL exported)
		endif (SVN_VERSION_OUTPUT MATCHES "Unversioned")
	endif (SVN_VERSION_RESULT)
endif (NOT GMTSAR_SOURCE_CODE_CONTROL_VERSION_STRING)

# The current GMTSAR version.
set (GMTSAR_VERSION_STRING "${GMTSAR_PACKAGE_NAME} ${GMTSAR_PACKAGE_VERSION_WITH_SVN_REVISION}")

set (GMTSAR_LONG_VERSION_STRING "${GMTSAR_PACKAGE_NAME} - ${GMTSAR_PACKAGE_DESCRIPTION_SUMMARY}, Version ${GMTSAR_PACKAGE_VERSION_WITH_SVN_REVISION}")

# use, i.e. don't skip the full RPATH for the build tree
set (CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# set the RPATH to be used when installing
if (NOT DEFINED GMTSAR_INSTALL_RELOCATABLE)
	set (GMTSAR_INSTALL_RELOCATABLE FALSE)
endif (NOT DEFINED GMTSAR_INSTALL_RELOCATABLE)
if (GMTSAR_INSTALL_RELOCATABLE)
	# make executables relocatable on supported platforms (relative RPATH)
	if (UNIX AND NOT CYGWIN)
		# find relative libdir from executable dir
		file (RELATIVE_PATH _rpath /${GMTSAR_BINDIR} /lib)
		# remove trailing /
		string (REGEX REPLACE "/$" "" _rpath "${_rpath}")
		if (APPLE)
			# relative RPATH on osx
			# CMP0042: CMake 3.0: MACOSX_RPATH is enabled by default
			set (CMAKE_MACOSX_RPATH ON)
			set (CMAKE_INSTALL_NAME_DIR @rpath)
			set (CMAKE_INSTALL_RPATH "@rpath;@executable_path/${_rpath}")
		else (APPLE)
			# relative RPATH on Linux, Solaris, etc.
			set (CMAKE_INSTALL_RPATH "\$ORIGIN/${_rpath}")
		endif (APPLE)
	endif (UNIX AND NOT CYGWIN)
else (GMTSAR_INSTALL_RELOCATABLE)
	# set absolute RPATH
	if (APPLE)
		# CMP0042: CMake 3.0: MACOSX_RPATH is enabled by default
		set (CMAKE_MACOSX_RPATH OFF)
		set (CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
	else (APPLE)
		# the RPATH to be used when installing, but only if it's not a
		# system directory
		list (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
			"${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
		if ("${isSystemDir}" STREQUAL "-1")
			set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
		endif ("${isSystemDir}" STREQUAL "-1")
	endif (APPLE)
endif (GMTSAR_INSTALL_RELOCATABLE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Make GNU and Intel C compiler default to C99
if (CMAKE_C_COMPILER_ID MATCHES "(GNU|Intel)" AND NOT CMAKE_C_FLAGS MATCHES "-std=")
	set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
endif ()

# vim: textwidth=78 noexpandtab tabstop=2 softtabstop=2 shiftwidth=2
