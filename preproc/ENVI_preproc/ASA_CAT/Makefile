include ../../../config.mk
PROG =  asa_cat
CSRCS = asa_cat.c 

OBJS =  $(CSRCS:.c=.o)
CLIBS = -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

clean:
	rm -f $(OBJS) tags core
