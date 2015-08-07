include ../../../config.mk
PROG =  ALOS_llt2rat
CSRCS = ALOS_llt2rat.c 

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
