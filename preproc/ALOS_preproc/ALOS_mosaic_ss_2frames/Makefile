include ../../../config.mk
PROG =  ALOS_mosaic_ss_2frames
CSRCS = ALOS_mosaic_ss_2frames.c 

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include -I../../../gmtsar $(GMT_INC)

CLIBS = -L../lib -lALOS -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(GMTSAR) $(GMT_LIB) $(INCLUDES)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(PROG) $(OBJS) tags core
