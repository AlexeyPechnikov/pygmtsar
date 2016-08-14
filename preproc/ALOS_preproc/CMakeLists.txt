#
#	$Id$
#

if (HAVE_M_LIBRARY)
	# link the math library
	set (LINK_LIBS m)
endif (HAVE_M_LIBRARY)

include_directories (include)

add_library (ALOS lib_src/utils.c lib_src/null_sio_struct.c
	lib_src/get_sio_struct.c lib_src/put_sio_struct.c lib_src/calc_dop.c
	lib_src/hermite_c.c lib_src/ALOS_ldr_orbit.c lib_src/ALOSE_orbits_utils.c
	lib_src/interpolate_ALOS_orbit.c lib_src/read_ALOS_sarleader.c
	lib_src/set_ALOS_defaults.c lib_src/write_ALOS_prm.c lib_src/rng_expand.c
	lib_src/rng_compress.c lib_src/rng_filter.c lib_src/find_fft_length.c
	lib_src/siocomplex.c lib_src/polyfit.c lib_src/plh2xyz.c lib_src/xyz2plh.c
	lib_src/cfft1d.c lib_src/swap32.c lib_src/swap16.c lib_src/fftpack.c
	include/image_sio.h include/lib_functions.h include/llt2xyz.h
	include/siocomplex.h)
target_link_libraries (ALOS ${LINK_LIBS})

set (LINK_LIBS ${LINK_LIBS} ALOS)

add_executable (ALOS_baseline ALOS_baseline/ALOS_baseline.c ALOS_baseline/ALOS_llt2rat_sub.c include/image_sio.h include/lib_functions.h include/llt2xyz.h)
target_link_libraries (ALOS_baseline ${LINK_LIBS})

add_executable (ALOS_fbd2fbs ALOS_fbd2fbs/ALOS_fbd2fbs.c include/image_sio.h include/lib_functions.h include/siocomplex.h)
target_link_libraries (ALOS_fbd2fbs ${LINK_LIBS})

add_executable (ALOS_fbd2fbs_SLC ALOS_fbd2fbs_SLC/ALOS_fbd2fbs_SLC.c include/image_sio.h include/lib_functions.h include/siocomplex.h)
target_link_libraries (ALOS_fbd2fbs_SLC ${LINK_LIBS})

add_executable (ALOS_fbd2ss ALOS_fbd2ss/ALOS_fbd2ss.c include/image_sio.h include/lib_functions.h include/siocomplex.h)
target_link_libraries (ALOS_fbd2ss ${LINK_LIBS})

add_executable (ALOS_fbs2fbd ALOS_fbs2fbd/ALOS_fbs2fbd.c include/image_sio.h include/siocomplex.h)
target_link_libraries (ALOS_fbs2fbd ${LINK_LIBS})

add_executable (ALOS_filter_SLC ALOS_filter_SLC/ALOS_filter_SLC.c include/image_sio.h include/siocomplex.h)
target_link_libraries (ALOS_filter_SLC ${LINK_LIBS})

add_executable (ALOS_llt2rat ALOS_llt2rat/ALOS_llt2rat.c include/image_sio.h include/lib_functions.h include/llt2xyz.h)
target_link_libraries (ALOS_llt2rat ${LINK_LIBS})

add_executable (ALOS_look ALOS_look/ALOS_look.c include/image_sio.h include/lib_functions.h include/llt2xyz.h)
target_link_libraries (ALOS_look ${LINK_LIBS})

add_executable (ALOS_merge ALOS_merge/ALOS_merge.c include/image_sio.h include/lib_functions.h)
target_link_libraries (ALOS_merge ${LINK_LIBS})

add_executable (ALOS_pre_process ALOS_pre_process/ALOS_pre_process.c ALOS_pre_process/parse_ALOS_commands.c ALOS_pre_process/read_ALOS_data.c ALOS_pre_process/read_ALOSE_data.c ALOS_pre_process/swap_ALOS_data_info.c ALOS_pre_process/roi_utils.c include/image_sio.h include/lib_functions.h)
target_link_libraries (ALOS_pre_process ${LINK_LIBS})

add_executable (ALOS_pre_process_SLC ALOS_pre_process_SLC/ALOS_pre_process_SLC.c ALOS_pre_process_SLC/parse_ALOS_commands.c ALOS_pre_process_SLC/read_ALOS_data_SLC.c ALOS_pre_process_SLC/swap_ALOS_data_info.c include/image_sio.h include/lib_functions.h)
target_link_libraries (ALOS_pre_process_SLC ${LINK_LIBS})

add_executable (ALOS_pre_process_SS ALOS_pre_process_SS/ALOS_pre_process_SS.c ALOS_pre_process_SS/parse_ALOS_commands.c ALOS_pre_process_SS/read_ALOS_data_SS.c ALOS_pre_process_SS/swap_ALOS_data_info.c include/image_sio.h include/lib_functions.h)
target_link_libraries (ALOS_pre_process_SS ${LINK_LIBS})

# add the install targets
install (TARGETS ALOS ALOS_baseline ALOS_fbd2fbs ALOS_fbd2fbs_SLC ALOS_fbd2ss ALOS_fbs2fbd ALOS_filter_SLC ALOS_llt2rat ALOS_look ALOS_merge ALOS_pre_process ALOS_pre_process_SLC ALOS_pre_process_SS
	ARCHIVE DESTINATION lib
	COMPONENT Runtime
	LIBRARY DESTINATION lib
	COMPONENT Runtime
	RUNTIME DESTINATION bin
	COMPONENT Runtime)

# vim: textwidth=78 noexpandtab tabstop=2 softtabstop=2 shiftwidth=2
