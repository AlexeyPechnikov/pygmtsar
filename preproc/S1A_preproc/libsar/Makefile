include ../../../config.mk
# Makefile for lib_src 
LIB = libsar.a
INCLUDES = -I../include -I$(BOOST_DIR)/include

SRCS =	doppler.cpp \
	gmtsar.cpp \
	instrument.cpp \
	orbit.cpp \
	planet.cpp \
	platform.cpp \
	scene.cpp \
	statevector.cpp \
	utilities.cpp

OBJS= $(SRCS:.cpp=.o)

all:	$(LIB)

install:	all

$(LIB) : $(OBJS)
	$(AR) rv $(LIB) $?
	$(RANLIB)  $(LIB)

clean:
	rm -f *.o $(LIB)

spotless:	clean
