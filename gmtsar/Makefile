#	$Id: Makefile 109 2015-01-19 23:01:24Z sandwell $
#
#	makefile for gmtsar directory (sbas.c not included yet)

include ../config.mk

INCLUDES	= $(GMT_INC) -Iinclude

PROGS_C		= bperp.c calc_dop_orb.c conv.c esarp.c offset_topo.c phase2topo.c \
		  phasediff.c phasefilt.c resamp.c xcorr.c extend_orbit.c \
		  SAT_llt2rat.c SAT_look.c SAT_baseline.c sbas.c make_gaussian_filter.c
PROGS_O         = $(PROGS_C:.c=.o)
PROGS           = $(PROGS_C:.c=)

LIB  		= lib_src/libgmtsar.$(LIBEXT)

lib_src/libgmtsar.$(LIBEXT):
		cd lib_src; $(MAKE) lib

#-------------------------------------------------------------------------------
#	software targets
#-------------------------------------------------------------------------------

all:		$(PROGS)

install:	all
		$(INSTALL) -d $(bindir)
		$(INSTALL) $(PROGS) $(bindir)
		ls csh/*.csh csh/*.sh > /tmp/t.lis
		$(INSTALL) `cat /tmp/t.lis` $(bindir)
		rm -f /tmp/t.lis

uninstall:
		cd $(bindir); rm -f $(PROGS)
		cd csh; ls *.csh > /tmp/t.lis
		cd $(bindir); rm -f `cat /tmp/t.lis`
		rm -f /tmp/t.lis

spotless::	clean
		cd lib_src; $(MAKE) spotless

clean:	
		rm -f *.o *% core tags $(PROGS)
		cd lib_src; $(MAKE) clean

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

$(PROGS):	$(PROGS_O) $(LIB)
		$(CC) $(LDFLAGS) $@.o $(GMTSAR) $(GMT_LIB) $(LALIBS) $(LIBS) -o $@
