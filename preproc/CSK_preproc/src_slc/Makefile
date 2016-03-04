include ../../../config.mk
PROG =  make_slc_csk
CSRCS = make_slc_csk.c

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include -I../../../gmtsar/include $(HDF5_CPPFLAGS)
CLIBS = -L../../../gmtsar/lib_src -lgmtsar -L../lib -lxmlC $(HDF5_LDFLAGS) $(HDF5_LIBS) -lm

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(CLIBS) 

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

clean:
	rm -f $(OBJS) tags core
