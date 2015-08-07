#       $Id$
#
#       makefile for top GMT5SAR directory

sinclude config.mk

# Currently, S1A must happen before the CSK, TSZ, and RS2 builds due to dependencies via links
# We will fix this so one can make anyting in any order

PREPROCESSORS	= ALOS ENVI ERS S1A CSK TSX RS2
DIRS		= gmtsar snaphu/src

all:	main preprocess

help::
		@grep '^#!' Makefile | cut -c3-
#!----------------- MAKE HELP FOR GMT5SAR-----------------
#!
#!make <target>, where <target> can be:
#!
#!all                : Compile everything
#!install            : Compile & install everything
#!clean              : Clean up and remove objects and library files
#!clean              : Clean and remove configured files as well
#!   ------------- Additional Targets   ------------- 
#!preprocess         : Only compile the preprocesors
#!main               : Only compile the main programs
#!install-preprocess : Only compile & install the preprocesors
#!install-main       : Only compile & install the main programs
#!

preprocess:
	for d in $(PREPROCESSORS); do \
		(cd preproc/$${d}_preproc; $(MAKE) all); \
	done

main:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) all); \
	done

install:	install-main install-preproc

install-preproc:
	for d in $(PREPROCESSORS); do \
		(cd preproc/$${d}_preproc; $(MAKE) install); \
	done
	$(INSTALL) preproc/ERS_preproc/scripts/virgin.PRM $(sharedir)
	$(INSTALL) preproc/ENVI_preproc/scripts/virgin_envisat.PRM $(sharedir)

install-main:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) install); \
	done
	$(INSTALL) -d $(sharedir)
	$(INSTALL) -d $(sharedir)/filters
	$(INSTALL) -d $(sharedir)/snaphu/config
	$(INSTALL) gmtsar/filters/[bfgsxy]* $(sharedir)/filters
	$(INSTALL) gmtsar/csh/snaphu.conf.* $(sharedir)/snaphu/config

uninstall:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) uninstall); \
	done
	rm -rf $(sharedir)

clean:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) clean); \
	done
	for d in $(PREPROCESSORS); do \
		(cd preproc/$${d}_preproc; $(MAKE) clean); \
	done

spotless:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) spotless); \
	done
	for d in $(PREPROCESSORS); do \
		(cd preproc/$${d}_preproc; $(MAKE) spotless); \
	done
	rm -rf $(sharedir) bin share
	rm -rf config.log config.status config.mk configure autom4te.cache
