include ../../../config.mk
PROG =  make_slc_csk
CSRCS = make_slc_csk.c

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include -I../../../gmtsar $(HDF5_CPPFLAGS)
CLIBS = -L../../../gmtsar -lgmtsar -L../lib -lxmlC $(HDF5_LDFLAGS) $(HDF5_LIBS) -lm

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(CLIBS) 

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
