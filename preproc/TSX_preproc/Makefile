include ../../config.mk
#
#   make all the files
#

DIRS	= libtinyxpath libsar src

all:
	ln -sf ../S1A_preproc/include include
	ln -sf ../S1A_preproc/libsar libsar
	ln -sf ../S1A_preproc/libtinyxpath libtinyxpath
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
