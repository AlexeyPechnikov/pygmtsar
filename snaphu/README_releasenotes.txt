Notable changes in v2.0.4 since v2.0.3:
---------------------------------------

* Fixed logic error that would sometimes cause crash in MCF solver
  when flow went from ground node to lower right corner node.

* Untabify source files (whitespace changes only)


Notable changes in v2.0.3 since v2.0.2:
---------------------------------------

* Fixed limit-checking error that would sometimes cause segmentation fault when
  width of masked area was only one pixel at right or bottom edges.
  

Notable changes in v2.0.2 since v2.0.1:
---------------------------------------

* Fixed indexing error that would sometimes cause segmentation fault when
  discharging boundary of masked region.
  

Notable changes in v2.0.1 since v2.0.0:
---------------------------------------

* Bug fixes affecting externally generated costs read from a file and
  handling of zero-cost primary and secondary arcs.
  

Notable changes in v2.0 since v1.4.2:
-------------------------------------

* The new -S option invokes behavior whereby snaphu will first run in
  tile mode to produce an unwrapped solution then feed this unwrapped
  solution back into the cost calculator and optimizer as if an
  unwrapped phase file were read from the input.  This option is
  equivalent to running snaphu in tile mode, then running snaphu again
  using the tile-mode output as an unwrapped input file using the -u
  option.  Tile parameters must be specified when using the -S option.

* Masking of input pixels is now supported.  A typical usage of
  masking would be in unwrapping repeat-pass interferograms where
  water areas are expected to have zero correlation; a water mask can
  be used to exclude such areas in order to reduce the execution time.
  
  When computing costs, arcs through masked areas are assumed to have
  zero cost, so nodes internal to masked areas can be ignored in the
  solver, thereby reducing the total number of nodes and usually
  decreasing the solution time.  An unwrapped phase output value will
  still be produced for each masked pixel, but the unwrapped value
  will not generally be reliable; often it will be whatever value is
  left from the initialization.

  The user can indicate masking by setting the magnitudes of the
  complex values of pixels to be masked to zero in the input file or
  by specifying a separate binary mask via the new -M option or the
  new BYTEMASKFILE keyword.  See the comments in the template
  configuration file for the BYTEMASKFILE keyword for a full
  description of the file format.  A fixed number of pixels at each
  edge of the input file can also be masked via the EDGEMASKTOP,
  EDGEMASKBOT, EDGEMASKLEFT, and EDGEMASKRIGHT keywords.

  If masked pixels separate regions of unmasked pixels, the unwrapped
  phase relationship between the disconnected regions will not
  generally be reliable, though a solution for each region will still
  produced (also see the NCONNNODEMIN keyword).

* The new -C command-line option takes a string argument that is
  parsed as a configuration line like one that would be put into a
  configuration file.  The string argument to the -C option may need
  to be protected from the shell by quotes to preserve whitespace.

* The integer types of several internal variables have changed now
  that 64-bit systems have generally replaced 32-bit systems.  The
  memory footprint of snaphu v2.0 is approximately 20% smaller than
  that of v1.4.2 when both are compiled as 64 bit binaries on an Intel
  system.  The memory footprint of snaphu v2.0 is similar to (about 2%
  larger than) that of v1.4.2 when both are compiled as a 32-bit
  binary, but v2.0 should be less constrained in terms of the input
  interferogram size.  The memory footprint of v2.0 when compiled as a
  64-bit binary is about 30% larger than when the same code (v2.0) is
  compiled as a 32-bit binary, although it is expected that most users
  will compile snaphu as a 64-bit binary nonetheless given the
  availability of memory on most systems.

* The syntax of the ASSEMBLEONLY keyword and the --assemble
  command-line option have changed so that they are now boolean flags.
  The name of the tile directory to assemble is specified through the
  new keyword TILEDIR or the --tiledir command-line option.
  Additionally, the new keyword DOTILEMASK allows the user to unwrap
  only selected tiles for manual experimentation and assembly of
  tiles.

* When unwrapping in tile mode, the default behavior is now to remove
  temporary tile files rather than to keep them.  Users can specify
  the new -k option on the command line or set the RMTMPTILE keyword
  to FALSE in a config file to keep temporary tile files.

* A warning is now displayed if running in tile mode and the tile
  overlap is less than a semi-arbitrary constant warning threshold.  A
  suggestion is displayed about increasing overlap or tile size if
  edge artifacts are present if the overlap is less than the tile size.

* Connected components are now supported in tile mode.  The connected
  components are assembled from the connected components of individual
  tiles and renumbered to be unique, but they will break at tile
  boundaries.  It is possible that a component will not be fully
  connected in the assembled output if the component was connected
  only in the tile overlap area that was discarded while assembling
  the full output.

* The connected components output file may be written as either 1-byte
  unsigned integers or 4-byte (32-bit) unsigned integers via the
  CONNCOMPOUTTYPE keyword.  Using 4-byte values allows more than 255
  connected components to be defined, although the default is still to
  use 1-byte values for backward compatibility.

* Support for Lp-norm cost functions has been added.  Note, however,
  that congruence is still enforced, meaning that the unwrapped phase
  will differ from the wrapped phase by only integer numbers cycles.
  Therefore, the use of an L2-norm cost functions is not equivalent to
  least-squares phase unwrapping approaches that do not enforce
  congruence (as described by Ghiglia and Romero, 1996).

* Experimental code for tree pruning has been included.  This behavior
  is controlled by the new NMAJORPRUNE and PRUNECOSTTHRESH keywords.
  This code has not been well tested and users are advised to use this
  code with only the lowest of expectations.

* The SOURCEMODE keyword and associated functionality have been removed
  due to changes in the source selection algorithm to allow the
  unwrapping of multiple regions that are disconnected by masked
  pixels.

* There have been many internal code cleanups that should be
  transparent to the user.

* Several miscellaneous bugs have been fixed and several minor
  enhancements have been added.
