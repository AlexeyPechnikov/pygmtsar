#       $Id: Makefile 236 2015-08-24 07:09:14Z pwessel $
#
#       makefile for top GMTSAR directory

sinclude config.mk

# Currently, S1A must happen before the CSK, TSZ, and RS2 builds due to dependencies via links
# We will fix this so one can make anyting in any order
#
PREPROCESSORS	= ALOS ERS S1A CSK TSX RS2 ENVI
#
DIRS		= gmtsar snaphu/src
ORBITS_URL	= http://topex.ucsd.edu/gmtsar/tar/ORBITS.tar
ORBITS		= ORBITS.tar

all:	main preprocess

help::
		@grep '^#!' Makefile | cut -c3-
#!----------------- MAKE HELP FOR GMTSAR-----------------
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
#!install-orbits     : Obtain ORBITS.tar, prompt for directory, and install files
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

install-orbits:
	wget $(ORBITS_URL) 2>/dev/null || curl -O $(ORBITS_URL)
	@if [ -f $(ORBITS) ]; then \
		read -p "==> Enter directory where the ORBITS subdirectory should be installed: " ORB_DIR ; \
		mkdir -p $$ORB_DIR/ORBITS ; \
		echo "Extracting orbits files into $$ORB_DIR/ORBITS";  \
		tar xf $(ORBITS) -C$$ORB_DIR/ORBITS ; \
		echo "Run configure --with-orbits-dir=$$ORB_DIR/ORBITS";  \
	else \
		echo "Unable to obtain $(ORBITS_URL) - Perhaps you have neither curl nor wget installed?";  \
	fi

uninstall:
	for d in $(DIRS) $(PREPROCESSORS); do \
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
