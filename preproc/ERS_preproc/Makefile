include ../../config.mk

DIRS	= ers_line_fixer read_data_file_dpaf read_data_file_ccrs read_sarleader_dpaf
SCRIPT	= scripts

all:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) all); \
	done

install:	all
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) install); \
	done
	$(INSTALL) $(SCRIPT)/dump_orbit_ers.pl $(bindir)
	$(INSTALL) $(SCRIPT)/ERS_pre_process $(bindir)
	$(INSTALL) $(bindir)/ENVI_baseline $(bindir)/ERS_baseline
	$(INSTALL) $(bindir)/ENVI_llt2rat  $(bindir)/ERS_llt2rat

uninstall:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) uninstall); \
	done
	rm -f $(bindir)/dump_orbit_ers.pl
	rm -f $(bindir)/ERS_pre_process
	rm -f $(bindir)/ERS_baseline
	rm -f $(bindir)/ERS_llt2rat

clean:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) clean); \
	done

spotless:	clean
