#-------------------------------------------------------------------------------
#  $Id: GNUmakefile 210 2015-08-07 15:46:10Z pwessel $
#
#       makefile for top GMT5SAR directory
#
#	!!! THIS MAKEFILE IS GMT5SAR GMT DEVELOPERS ONLY !!!
#
#	This makefile extends the regular Makefile by adding
#	targets to build tarballs, get latest GNU config files.
#
#	Date:		27-JAN-2013
#-------------------------------------------------------------------------------

include Makefile

# Where to place the gmtsar tarball
#FTPDIR	= /topex/ftp/pub/gmtsar
FTPDIR	= ..

create:
# We make the $(PACKAGE_TARNAME)-$(GMT5SAR_VERSION) link from scratch each time
	cd ..; rm -f $(PACKAGE_TARNAME)-$(GMT5SAR_VERSION); $(LN_S) $(notdir $(PWD)) $(PACKAGE_TARNAME)-$(GMT5SAR_VERSION)

tar:	distro

distro: create
	$(MAKE) spotless
	find . -type l -exec rm -f {} \;
	COPYFILE_DISABLE=true $(GNUTAR) --owner 0 --group 0 --mode a=rX,u=rwX -cjhvf $(FTPDIR)/$(PACKAGE_TARNAME)-$(GMT5SAR_VERSION)-src.tar.bz2 \
		-C .. $(PACKAGE_TARNAME)-$(GMT5SAR_VERSION) --exclude GNUmakefile --exclude .svn
	rm -f $(PACKAGE_TARNAME)-$(GMT5SAR_VERSION)

latest-config:
		curl "http://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub;hb=HEAD" -s -R -o config.sub
		curl "http://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess;hb=HEAD" -s -R -o config.guess
		curl "http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_blas.m4" -s -R -o ax_blas.m4
		curl "http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_lapack.m4" -s -R -o ax_lapack.m4
