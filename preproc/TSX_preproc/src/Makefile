include ../../../config.mk
PROG =  make_slc_tsx
CSRCS = make_slc_tsx.cpp tsx.cpp Burst.cpp Cosar.cpp Header.cpp

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

clean:
	rm -f $(OBJS) tags core
