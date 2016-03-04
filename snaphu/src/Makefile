#	$Id: Makefile 33 2013-04-06 05:37:15Z pwessel $
#
#	makefile for snaphu directory
#
# Makefile for statistical-cost network-flow algorithm for phase unwrapping 
# Curtis W. Chen
# Copyright 2002 Board of Trustees, Leland Stanford Jr. University
#
# If you want to copy the man page and executable to a system directory, 
# specify the appropriate directories in the INSTALLDIR and MANDIR 
# variables and do 'make install' (you will probably need to be root).
#
# If you specify -D NO_CS2, the program will be compiled without the
# CS2 MCF solver module.

include ../../config.mk

PROGS_C		= snaphu.c
LIB_C		= snaphu_tile.c \
		  snaphu_solver.c \
		  snaphu_io.c \
		  snaphu_util.c \
		  snaphu_cost.c \
		  snaphu_cs2.c 
PROGS_O         = $(PROGS_C:.c=.o)
PROGS           = $(PROGS_C:.c=)
MANDIR		= /usr/local/man
SNAPHUMAN	= ../man/man1/snaphu.1
LIB_O         = $(LIB_C:.c=.o)

#-------------------------------------------------------------------------------
#	software targets
#-------------------------------------------------------------------------------

all:		$(PROGS)

install:	all
		$(INSTALL) -d $(bindir)
		$(INSTALL) $(PROGS) $(bindir)

install-man:
		cp $(SNAPHUMAN) $(MANDIR)/man1

uninstall:
		cd $(bindir); rm -f $(PROGS)

spotless::	clean

clean:	
		rm -f *.o *% core tags $(PROGS)

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

$(PROGS):	$(PROGS_O) $(LIB_O)
		$(CC) $(LDFLAGS) $@.o $(LIB_O) $(LIBS) -o $@
