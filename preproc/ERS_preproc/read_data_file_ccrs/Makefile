include ../../../config.mk
PROG =  read_data_file_ccrs
CSRCS = read_data_file_ccrs.c

OBJS =  $(CSRCS:.c=.o)	

INCLUDES = -I../include
LIBS = -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(LIBS) 

all:	$(PROG)

clean:
	rm -f $(OBJS) tags core *.c.*

install:
	$(INSTALL) $(PROG) $(bindir)/$(PROG)

uninstall:
	rm -f $(bindir)/$(PROG)
