include ../../../config.mk
PROG =  make_raw_csk
CSRCS = make_raw_csk.c

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include -I../../../gmtsar $(HDF5_CPPFLAGS)
CLIBS = -L../../../gmtsar -lgmtsar -L../lib -lxmlC $(HDF5_LDFLAGS) $(HDF5_LIBS) -lm

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(CLIBS) 

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
