include ../../config.mk

LIB	= lib/ibENVI.a
DIRS	= ASA_CAT Dop_orbit ENVI_baseline ENVI_decode ENVI_llt2rat ENVI_look
SCRIPT	= scripts

all:	$(LIB)
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) all); \
	done

install:	all
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) install); \
	done
	$(INSTALL) $(SCRIPT)/dump_orbit_envi.pl $(bindir)
	$(INSTALL) $(SCRIPT)/dump_time_envi.pl $(bindir)
	$(INSTALL) $(SCRIPT)/find_auxi.pl $(bindir)
	$(INSTALL) $(SCRIPT)/ENVI_pre_process $(bindir)

uninstall:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) uninstall); \
	done
	rm -f $(bindir)/dump_orbit_envi.pl
	rm -f $(bindir)/dump_time_envi.pl
	rm -f $(bindir)/find_auxi.pl
	rm -f $(bindir)/ENVI_pre_process

clean:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) clean); \
	done
	(cd lib_src; $(MAKE) clean)

$(LIB):
	mkdir -p lib
	(cd lib_src; $(MAKE) all)

spotless:	clean
	rm -rf lib
