include ../../../config.mk
PROG =  ENVI_baseline
CSRCS = ENVI_baseline.c \
	ENVI_llt2rat_sub.c

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../../ALOS_preproc/include
CLIBS = -L../../ALOS_preproc/lib -lALOS -L../lib -lENVI -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

clean:
	rm -f $(OBJS) tags core
