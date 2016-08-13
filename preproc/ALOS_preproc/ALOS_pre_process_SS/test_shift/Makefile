include ../../../config.mk
PROG =  test_shift
CSRCS = test_shift.c \
	shift.c
OS = `arch`

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../../include
CLIBS = -L../../lib/$(OS) -lALOS -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

install:
	$(INSTALL) $(PROG) ../bin/$(OS)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
