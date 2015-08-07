include ../../../config.mk
PROG =  ALOS_look
CSRCS = ALOS_look.c 

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include

CLIBS = -L../lib -lALOS -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

clean:
	rm -f $(OBJS) tags core
