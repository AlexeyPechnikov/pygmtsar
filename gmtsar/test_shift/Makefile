PROG =  test_shift
CSRCS = test_shift.c 

FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC =  /opt/local/include/gmt 

CLIBS = -L/opt/local/lib -lgmt -L../lib_src -lgmtsar -lm
FLIBS = 
CFLAGS = -O2 -Bstatic $(M64) -w

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS) $(VEC)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC) -I/usr/local/GMT5SAR/gmtsar/include

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
