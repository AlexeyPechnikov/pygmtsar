#	$Id: Makefile 64 2013-04-16 04:41:43Z pwessel $
#
#	makefile for gmtsar/lib_src directory
#
include ../../config.mk

INCLUDES	= $(GMT_INC) -I../include

LIB_C		= aastretch.c acpatch.c conv2d.c do_freq_xcorr.c \
		  do_time_int_xcorr.c fft_bins.c fft_interpolate_routines.c \
		  file_stuff.c geoxyz.c get_locations.c get_params.c \
		  highres_corr.c intp_coef.c parse_xcorr_input.c plxyz.c \
		  print_results.c radopp.c read_xcorr_data.c rmpatch.c \
		  rng_cmp.c rng_ref.c shift.c sio_struct.c siocomplex.c \
		  spline.c trans_col.c utils.c utils_complex.c set_prm_defaults.c\
		  ldr_orbit.c read_orb.c interpolate_orbit.c polyfit.c calc_dop.c

LIB_O		= $(LIB_C:.c=.o)
LIB		= libgmtsar.$(LIBEXT)

#-------------------------------------------------------------------------------
#	library
#-------------------------------------------------------------------------------

$(LIB):		$(LIB_O)
		$(AR) cvur $@ $?
		$(RANLIB) $@

lib:		$(LIB)
		$(INSTALL) -d ../lib
		$(INSTALL) -m 644 $(LIB:.$(LIBEXT)=.a) ../lib
		if [ ! $(LIBEXT) = a ]; then $(INSTALL) -m 755 $(LIB) ../lib; fi

clean:
		rm -f *% *.o *.a

spotless::	clean
		rm -rf ../lib
