#
#	$Id$
#

add_executable (snaphu
	snaphu.c
	snaphu_cost.c
	snaphu_cs2.c
	snaphu_io.c
	snaphu_solver.c
	snaphu_tile.c
	snaphu_util.c
	snaphu_cs2types.h
	snaphu.h)

if (HAVE_M_LIBRARY)
	# link the math library
	target_link_libraries (snaphu m)
endif (HAVE_M_LIBRARY)

# add the install targets
install (TARGETS snaphu
	ARCHIVE DESTINATION lib
	COMPONENT Runtime
	LIBRARY DESTINATION lib
	COMPONENT Runtime
	RUNTIME DESTINATION bin
	COMPONENT Runtime)

# vim: textwidth=78 noexpandtab tabstop=2 softtabstop=2 shiftwidth=2
