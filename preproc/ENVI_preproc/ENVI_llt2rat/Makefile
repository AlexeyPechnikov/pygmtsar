include ../../../config.mk
PROG =  ENVI_llt2rat
CSRCS = ENVI_llt2rat.c 

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../../ALOS_preproc/include

CLIBS = -L../../ALOS_preproc/lib -lALOS -L../lib -lENVI -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
