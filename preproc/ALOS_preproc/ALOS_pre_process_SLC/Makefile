include ../../../config.mk
PROG =  ALOS_pre_process_SLC
CSRCS = ALOS_pre_process_SLC.c \
	parse_ALOS_commands.c \
	read_ALOS_data_SLC.c \
	swap_ALOS_data_info.c 


OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include
CLIBS = -L../lib -lALOS -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
