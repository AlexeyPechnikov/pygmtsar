#
#	$Id: CMakeLists.txt 15759 2016-02-24 23:36:11Z pwessel $
#
# To modify the cmake process: Edit your cmake/ConfigUser.cmake file
#
# To build out-of-source do (example):
#
#	mkdir build
#	cd build
#	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
#
# CMAKE_BUILD_TYPE can be: empty, Debug, Release, RelWithDebInfo or MinSizeRel
#

# Make sure the user doesn't play dirty with symlinks
get_filename_component (srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component (bindir "${CMAKE_BINARY_DIR}" REALPATH)

# Disallow in-source builds
if (${srcdir} STREQUAL ${bindir})
	message(FATAL_ERROR "In-source builds are not allowed. "
	"Please create a directory and run cmake from there, passing the path "
	"to this source directory as the last argument. This process created "
	"the file `CMakeCache.txt' and the directory `CMakeFiles' in ${srcdir}. "
	"Please remove them.")
endif (${srcdir} STREQUAL ${bindir})

# Define minimum CMake version required
cmake_minimum_required (VERSION 2.8.5)

# Use NEW behavior with newer CMake releases
foreach(p
		CMP0025 # CMake 3.0: Compiler id for Apple Clang is now AppleClang
		)
	if(POLICY ${p})
		cmake_policy(SET ${p} NEW)
	endif()
endforeach()

# Define project name and language
project (GMTSAR C)

# Where to find our CMake modules (this variable is visible in subdirectories).
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/"
	CACHE INTERNAL "Location of our custom CMake modules." FORCE)

# Include configuration options (default options and options overridden by user).
include (ConfigCMake)

# Find UNIX commands
include (FindUnixCommands)
find_program (SVN svn)
find_program (GS gs gswin64)

# Find dependent libraries
find_package (GMT REQUIRED)
find_package (TIFF REQUIRED)
find_package (LAPACK)

# Add subdirectories
add_subdirectory (gmtsar)

# Configuration done
message(
	"*  Options:\n"
	"*  GMT library                : ${GMT_LIBRARY}\n"
	"*  GMT include dir            : ${GMT_INCLUDE_DIR}\n"
	"*  TIFF library               : ${TIFF_LIBRARY}\n"
	"*  LAPACK library             : ${LAPACK_LIBRARIES} ${LAPACK_lapack_LIBRARY}\n"
	"*\n"
	"*  Locations:\n"
	"*  Installing GMTSAR in       : ${CMAKE_INSTALL_PREFIX}\n"
)

# For debugging: print all set variables
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#	message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

# vim: textwidth=78 noexpandtab tabstop=2 softtabstop=2 shiftwidth=2
