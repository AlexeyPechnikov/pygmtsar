include ../../../config.mk

INCLUDES = -I../include -I$(GMTSARHOME)/gmtsar -I$(TIFF_INC) $(GMT_INC)
CLIBS = -L../lib -lxmlC $(GMTSAR) $(GMT_LIB) $(LALIBS) -L$(TIFF_LIB) -ltiff -lm

assemble_tops: assemble_tops.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ assemble_tops.o $(CLIBS) 

all:	assemble_tops

install:
	$(INSTALL) assemble_tops $(bindir)

clean:
	rm -f *.o assemble_tops tags core
