include ../../../config.mk
PROG =  ALOS_fbd2ss
CSRCS = ALOS_fbd2ss.c

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
