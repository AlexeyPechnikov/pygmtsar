////
$Id: CMakeLists.txt 15759 2016-02-24 23:36:11Z pwessel $
////

Information for building GMTSAR with CMake
------------------------------------------

Quick start:
~~~~~~~~~~~~

Checkout GMTSAR from its subversion repository:

  $ svn co svn://gmtserver.soest.hawaii.edu/GMTSAR/trunk gmtsar-dev
  $ cd gmtsar-dev
  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
  $ make -j
  $ make -j install

Invoking CMake:
~~~~~~~~~~~~~~~

CMake supports out-of-source builds and puts generated files in a completely
separate directory, so that the source tree is unchanged. Trying to run CMake
in the source tree will result in an error. CMake will first check for the
minimal required version and then search for required and optional libraries
(GMT, TIFF, LAPACK) in the default locations. The configuration step will fail
if dependencies are not met.

Updating:
~~~~~~~~~

Assuming you did not delete the build directory and that your current
working directory is the build directory this is just as simple as

  $ cd ..
  $ svn up
  $ cd -
  $ make -jx install

CMake will detect any changes to the source files and will automatically
reconfigure. If you deleted all files inside the build directory you have to
run cmake again manually.

////
vim: ft=asciidoc tw=78 et ts=2 sts=2 sw=2 autoindent smartindent
////
