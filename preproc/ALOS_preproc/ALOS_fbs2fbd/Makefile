include ../../../config.mk
PROG =  ALOS_fbs2fbd
CSRCS = ALOS_fbs2fbd.c

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include

CLIBS = -L../lib -lALOS -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
