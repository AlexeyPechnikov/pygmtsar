include ../../../config.mk

INCLUDES = -I$(GMTSARHOME)/gmtsar/include $(GMT_INC)
CLIBS = $(GMTSAR) $(GMT_LIB) $(LALIBS) -lm

stitch_tops: stitch_tops.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ stitch_tops.o $(CLIBS) 

all:	stitch_tops

install:
	$(INSTALL) stitch_tops $(bindir)

clean:
	rm -f *.o stitch_tops tags core
