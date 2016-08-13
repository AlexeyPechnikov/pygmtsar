include ../../../config.mk
PROG =  read_sarleader_dpaf
CSRCS = read_sarleader_dpaf.c make_prm_dpaf.c

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
