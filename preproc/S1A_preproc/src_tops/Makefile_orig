include ../../../config.mk
PROG =  make_s1a_tops
CSRCS = make_s1a_tops.c

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include -I$(GMTSARHOME)/gmtsar/include -I$(TIFF_INC) $(GMT_INC)
CLIBS = -L../lib -lxmlC $(GMTSAR) $(GMT_LIB) $(LALIBS) -L$(TIFF_LIB) -ltiff -lm

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(CLIBS) 

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

clean:
	rm -f $(OBJS) tags core
