include ../../../config.mk
PROG =  ENVI_look
CSRCS = ENVI_look.c 

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INCLUDES = -I../../ALOS_preproc/include 

CLIBS =	-L../../ALOS_preproc/lib -lALOS -L../lib -lENVI -lm 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
