include ../../../config.mk
PROG =  asa_im_decode
CSRCS = asa_im_decode.c 
OBJS =  $(CSRCS:.c=.o)
CLIBS = -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
