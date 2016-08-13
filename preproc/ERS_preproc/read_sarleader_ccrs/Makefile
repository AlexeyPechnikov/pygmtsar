include ../../../config.mk
PROG =  read_sarleader_ccrs
CSRCS = read_sarleader_ccrs.c make_prm_ccrs.c

OBJS =  $(CSRCS:.c=.o)	

LIBS = -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(LIBS)

all:	$(PROG)

clean:
	rm -f $(OBJS) tags core

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)
