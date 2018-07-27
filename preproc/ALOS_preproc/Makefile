include ../../config.mk

DIRS	= ALOS_pre_process ALOS_merge ALOS_fbd2fbs \
	  ALOS_fbd2fbs_SLC ALOS_fbd2ss ALOS_baseline \
	  ALOS_pre_process_SLC ALOS_mosaic_ss ALOS_look \
	  ALOS_mosaic_ss ALOS_mosaic_ss_2frames ALOS_pre_process_SS
 
LIB	= lib/libALOS.a

all:	$(LIB)
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) all); \
	done

install:	all
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) install); \
	done

uninstall:
	for d in $(DIRS); do \
		(cd $$d; $(MAKE) uninstall); \
	done

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
