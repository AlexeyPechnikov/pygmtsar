include ../../../config.mk
# Makefile for lib_src 
LIB = libtinyxpath.a
INCLUDES = -I../include

SRCS =	action_store.cpp \
	lex_util.cpp \
	node_set.cpp \
	tinystr.cpp \
	tinyxml.cpp \
	tinyxmlerror.cpp \
	tinyxmlparser.cpp \
	tokenlist.cpp \
	xml_util.cpp \
	xpath_expression.cpp \
	xpath_processor.cpp \
	xpath_stack.cpp \
	xpath_static.cpp \
	xpath_stream.cpp \
	xpath_syntax.cpp \

OBJS= $(SRCS:.cpp=.o)

all:	$(LIB)

install:	all

$(LIB) : $(OBJS)
	$(AR) rv $(LIB) $?
	$(RANLIB) $(LIB)

clean:
	rm -f *.o $(LIB)

spotless:	clean
