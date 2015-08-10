include ../../../config.mk
PROG =  make_raw_csk
CSRCS = make_raw_csk.cpp csk.cpp

OBJS =  $(CSRCS:.cpp=.o)
INCLUDES = -I../include -I../libsar/include -I$(BOOST_DIR)/include $(HDF5_CPPFLAGS)
CLIBS = -L../libtinyxpath -ltinyxpath -L../libsar -lsar -L$(BOOST_DIR)/lib \
	-lboost_system$(MT) -lboost_program_options$(MT) \
	-lboost_filesystem$(MT) -lboost_date_time$(MT) $(HDF5_LDFLAGS) $(HDF5_LIBS)

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(CLIBS) 

make_raw_csk.o:
	$(CXX) $(CXXFLAGS) -c $(CSRCS) $(INCLUDES)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

clean:
	rm -f $(OBJS) tags core
