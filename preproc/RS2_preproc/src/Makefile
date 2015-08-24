include ../../../config.mk
PROG =  make_slc_rs2
CSRCS = make_slc_rs2.cpp rs2.cpp 

OBJS =  $(CSRCS:.cpp=.o)
INCLUDES = -I../include -I../libsar/include -I$(BOOST_INC) -I$(TIFF_INC)
CLIBS = -L../libtinyxpath -ltinyxpath -L../libsar -lsar -L$(BOOST_LIB) \
	-lboost_system$(MT) -lboost_program_options$(MT) \
	-lboost_filesystem$(MT) -lboost_date_time$(MT) -L$(TIFF_LIB) -ltiff -lm

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(CLIBS) 

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
