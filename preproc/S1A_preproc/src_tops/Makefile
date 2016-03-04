include ../../../config.mk

INCLUDES = -I../include -I$(GMTSARHOME)/gmtsar/include -I$(TIFF_INC) $(GMT_INC)
CLIBS = -L../lib -lxmlC $(GMTSAR) $(GMT_LIB) $(LALIBS) -L$(TIFF_LIB) -ltiff -lm

make_s1a_tops: make_s1a_tops.o
	$(CC) $(CFLAGS) -o $@ make_s1a_tops.o $(CLIBS) 

make_s1a_tops_6par: make_s1a_tops_6par.o
	$(CC) $(CFLAGS) -o $@ make_s1a_tops_6par.o $(CLIBS) 

all:	make_s1a_tops make_s1a_tops_6par

install:
	$(INSTALL) make_s1a_tops $(bindir)
	$(INSTALL) make_s1a_tops_6par $(bindir)

clean:
	rm -f *.o make_s1a_tops make_s1a_tops_6par tags core
