include ../../../config.mk

INCLUDES = -I../include -I$(GMTSARHOME)/gmtsar $(GMT_INC)
CLIBS = -L../lib -lxmlC $(GMTSAR) $(GMT_LIB) $(LALIBS) -lm

spectral_diversity: spectral_diversity.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ spectral_diversity.o $(CLIBS) 

all:	spectral_diversity

install:
	$(INSTALL) spectral_diversity $(bindir)

clean:
	rm -f *.o spectral_diversity tags core
