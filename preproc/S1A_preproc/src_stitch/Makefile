include ../../../config.mk

INCLUDES = -I$(GMTSARHOME)/gmtsar $(GMT_INC)
CLIBS = $(GMTSAR) $(GMT_LIB) $(LALIBS) -lm

stitch_tops: stitch_tops.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ stitch_tops.o $(CLIBS) 

merge_swath: merge_swath.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ merge_swath.o $(CLIBS)

all:	stitch_tops merge_swath

install:
	$(INSTALL) stitch_tops $(bindir)
	$(INSTALL) merge_swath $(bindir)

clean:
	rm -f *.o stitch_tops merge_swath tags core
