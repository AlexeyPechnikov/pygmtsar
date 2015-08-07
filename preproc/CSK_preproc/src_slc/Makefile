include ../../../config.mk
PROG =  make_slc_csk
CSRCS = make_slc_csk.cpp csk.cpp

OBJS =  $(CSRCS:.cpp=.o)
INCLUDES = -I../include -I../libsar/include -I$(BOOST_DIR)/include
CLIBS = -L../libtinyxpath -ltinyxpath -L../libsar -lsar -L$(BOOST_DIR)/lib \
	-lboost_system$(MT) -lboost_program_options$(MT) \
	-lboost_filesystem$(MT) -lboost_date_time$(MT) -lhdf5_cpp -lhdf5 -lm

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(CLIBS) 

make_slc_csk.o:
	$(CXX) $(CXXFLAGS) -c $(CSRCS) $(INCLUDES)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

clean:
	rm -f $(OBJS) tags core
