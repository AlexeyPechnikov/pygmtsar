### makefile for bert to create a lib of object code
#%// BK 13-Jun-2003
#%// MA 24-Dec-2008    -m32 switch for 64-bit machines
#%// MA 05-Jan-2010    64-bit native with new epr_api_2.2

SRC = 		epr_api.c \
		epr_band.c \
		epr_bitmask.c \
		epr_core.c \
		epr_dataset.c \
		epr_dddb.c \
		epr_dsd.c \
		epr_dump.c \
		epr_field.c \
		epr_msph.c \
		epr_param.c \
		epr_product.c \
		epr_ptrarray.c \
		epr_record.c \
		epr_string.c \
		epr_swap.c \
		epr_typconv.c

OBJ = 		$(SRC:.c=.o)
CC =		gcc
#CFLAGS = 	-O3
#CFLAGS = 	-m32         # for 64-bit systems, if epr_api version < 2.2

#CFLAGS = 	-fPIC -ansi -c

.c.o:	
		$(CC) $(CFLAGS) -c -o $(@) $<

default:	$(OBJ)
		echo "objects created done"
#		ar ...

clean:		
		rm *.o
		@echo "removed all object files"


