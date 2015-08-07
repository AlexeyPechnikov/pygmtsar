include ../../config.mk
#
#   make all the files
#
DIRS	= libtinyxpath libsar src

all:
	$(LN_S) -f ../S1A_preproc/include include
	$(LN_S) -f ../S1A_preproc/libsar libsar
	$(LN_S) -f ../S1A_preproc/libtinyxpath libtinyxpath
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) all); \
	done

install:	all
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) install); \
	done

clean:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) clean); \
	done

spotless:	clean
