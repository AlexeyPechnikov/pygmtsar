include ../../../config.mk
# Makefile for lib
LIB = libxmlC.a
INCLUDES = -I../include 

SRCS =	xml.c 

OBJS= $(SRCS:.c=.o)

all:	$(LIB)

install:	all

$(LIB) : $(OBJS)
	$(AR) rv $(LIB) $?
	$(RANLIB)  $(LIB)

clean:
	rm -f *.o $(LIB)

spotless:	clean
