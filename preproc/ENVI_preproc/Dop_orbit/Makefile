include ../../../config.mk
PROG =  calc_dop_orb_envi
CSRCS = calc_dop_orb_envi.c 

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../../ALOS_preproc/include

CLIBS = -L../../ALOS_preproc/lib -L../lib -lENVI -lALOS -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
