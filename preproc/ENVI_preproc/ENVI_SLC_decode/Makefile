include ../../../config.mk
PROG =  envi_slc_decode
CSRCS = envi_slc_decode.c


OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include -I../../../gmtsar/include -I../../../gmtsar -I../../S1A_preproc/include -I$(TIFF_INC)
CLIBS = -L../../../gmtsar -lgmtsar -L../lib -L../../S1A_preproc/lib  -lxmlC -lm

LIB_SRC =       epr_api-2.3/src/epr_api.c \
                epr_api-2.3/src/epr_band.c \
                epr_api-2.3/src/epr_bitmask.c \
                epr_api-2.3/src/epr_core.c \
                epr_api-2.3/src/epr_dataset.c \
                epr_api-2.3/src/epr_dddb.c \
                epr_api-2.3/src/epr_dsd.c \
                epr_api-2.3/src/epr_dump.c \
                epr_api-2.3/src/epr_field.c \
                epr_api-2.3/src/epr_msph.c \
                epr_api-2.3/src/epr_param.c \
                epr_api-2.3/src/epr_product.c \
                epr_api-2.3/src/epr_ptrarray.c \
                epr_api-2.3/src/epr_record.c \
                epr_api-2.3/src/epr_string.c \
                epr_api-2.3/src/epr_swap.c \
                epr_api-2.3/src/epr_typconv.c


all:	BEST envi_slc_decode envisat_dump_data envisat_dump_header

### create object code BEST library functions ###
BEST:		$(LIB_SRC)
	(cd epr_api-2.3/src/; make)

install:
	$(INSTALL) $(PROG) $(bindir)
	$(INSTALL) envisat_dump_data $(bindir)
	$(INSTALL) envisat_dump_header $(bindir)

clean:
	rm -f $(OBJS) envi_slc_decode envisat_dump_data envisat_dump_header tags core
	(cd epr_api-2.3/src/; make clean)

envi_slc_decode:
		$(CC) $(CFLAGS) $(INCLUDES) $(CLIBS)$(LFLAGS) \
		-o envi_slc_decode \
		epr_api-2.3/src/*.o \
		envi_slc_decode.c

envisat_dump_data:
		$(CC) $(CFLAGS) $(INCLUDES) $(CLIBS)$(LFLAGS) \
		-o envisat_dump_data \
		epr_api-2.3/src/*.o \
		envisat_dump_data.c

envisat_dump_header:
		$(CC) $(CFLAGS) $(INCLUDES) $(CLIBS)$(LFLAGS) \
		-o envisat_dump_header \
		epr_api-2.3/src/*.o \
		envisat_dump_header.c
