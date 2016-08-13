include ../../config.mk
#
#   make all the files
#
DIRS	= src_slc src_raw

all:
	$(LN_S) -f ../S1A_preproc/include include
	$(LN_S) -f ../S1A_preproc/lib lib
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) all); \
	done

install:	all
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) install); \
	done

uninstall:
	for d in src_slc src_raw; do \
		(cd $$d; $(MAKE) uninstall); \
	done

clean:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) clean); \
	done

spotless:	clean
