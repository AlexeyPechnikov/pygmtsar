include ../../config.mk
#
#   make all the files
#
DIRS	= lib src_swath src_tops src_orbit src_stitch src_spec_div

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
	(cd src; $(MAKE) uninstall)

clean:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) clean); \
	done

spotless:	clean
