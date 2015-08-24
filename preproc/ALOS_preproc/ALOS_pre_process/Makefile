include ../../../config.mk
PROG =  ALOS_pre_process
CSRCS = ALOS_pre_process.c \
	parse_ALOS_commands.c \
	read_ALOS_data.c \
	read_ALOSE_data.c \
	swap_ALOS_data_info.c \
	roi_utils.c

OBJS =  $(CSRCS:.c=.o)
INCLUDES = -I../include
CLIBS = -L../lib -lALOS -lm

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)  $(FLIBS)

all:	$(PROG)

install:
	$(INSTALL) $(PROG) $(bindir)

uninstall:
	rm -f $(bindir)/$(PROG)

clean:
	rm -f $(OBJS) tags core
