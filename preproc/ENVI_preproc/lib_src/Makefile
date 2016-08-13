include ../../../config.mk
# Makefile for lib src 
LIB = libENVI.a
INCLUDES = -I../../ALOS_preproc/include

SRCS =	read_ENVI_orb.c \
	ENVI_ldr_orbit.c

OBJS= $(SRCS:.c=.o)

$(LIB) : $(OBJS)
	$(AR) r $(LIB) $?
	$(RANLIB) $(LIB)

all:	install

install:	$(LIB)
	$(INSTALL) *.a ../lib

clean:
	rm -f *.a *.o 
