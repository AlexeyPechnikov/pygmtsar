sinclude ../../../config.mk
# Makefile for lib src 
LIB = libALOS.a
INCLUDES = -I../include

SRCS =	utils.c \
	null_sio_struct.c \
	get_sio_struct.c \
	put_sio_struct.c \
	calc_dop.c \
	hermite_c.c \
	ALOS_ldr_orbit.c \
	ALOSE_orbits_utils.c \
	interpolate_ALOS_orbit.c \
	read_ALOS_sarleader.c \
	set_ALOS_defaults.c \
	write_ALOS_prm.c \
	rng_expand.c \
	rng_compress.c \
	rng_filter.c \
	find_fft_length.c \
	siocomplex.c \
	polyfit.c \
	plh2xyz.c \
	xyz2plh.c \
	cfft1d.c \
	swap32.c \
	swap16.c \
	fftpack.c

OBJS= $(SRCS:.c=.o)

$(LIB) : $(OBJS)
	$(AR) r $(LIB) $?
	$(RANLIB) $(LIB)

all:	install

install:	$(LIB)
	$(INSTALL) *.a ../lib

clean:
	rm -f *.a *.o 
