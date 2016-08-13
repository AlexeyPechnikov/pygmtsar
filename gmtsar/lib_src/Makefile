#	$Id: Makefile 109 2015-01-19 23:01:24Z sandwell $
#
#    makefile for gmtsar directory
#
include ../../config.mk

INCLUDES	= $(GMT_INC) -I../include

LIB_C		= aastretch.c acpatch.c calc_dop.c conv2d.c do_freq_xcorr.c \
		  do_time_int_xcorr.c fft_bins.c fft_interpolate_routines.c \
		  file_stuff.c geoxyz.c get_locations.c get_params.c hermite_c.c \
		  highres_corr.c interpolate_orbit.c intp_coef.c ldr_orbit.c \
		  parse_xcorr_input.c plxyz.c polyfit.c print_results.c radopp.c \
		  read_orb.c read_xcorr_data.c SAT_llt2rat_sub.c \
		  rmpatch.c rng_cmp.c rng_ref.c set_prm_defaults.c shift.c \
		  sio_struct.c siocomplex.c spline.c trans_col.c utils.c utils_complex.c \
		  write_orb.c sbas_utils.c

LIB_O		= $(LIB_C:.c=.o)
LIB		= libgmtsar.$(LIBEXT)

#-------------------------------------------------------------------------------
#	library
#-------------------------------------------------------------------------------

$(LIB):		$(LIB_O)
		$(AR) cvur $@ $?
		$(RANLIB) $@

lib:		$(LIB)

clean:
		rm -f *% *.o *.a

spotless::	clean
