#!/bin/bash
#	$Id$
# Uninstall GMT5SAR installation (share and executables)
#
printf "gmtsar_uninstall.sh will uninstall the GMT5SAR executables and share directory\n" >&2
printf "You may need sudo privileges on this computer.\n\nContinue? (y/n) [y]:" >&2
read answer
if [ "X$answer" = "Xn" ]; then
        exit 0
fi
here=`pwd`
sharedir=`gmtsar_sharedir.csh`
bin=`which gmtsar_sharedir.csh`
dir=`dirname $bin`
modules="ALOS_fbd2fbs ALOS_fbd2fbs_SLC ALOS_fbd2ss ALOS_look ALOS_merge \
	ALOS_pre_process ALOS_pre_process_SLC ALOS_pre_process_SS ENVI_baseline ENVI_llt2rat ENVI_look \
	ENVI_pre_process ERS_baseline ERS_llt2rat ERS_pre_process SAT_baseline SAT_llt2rat SAT_look \
	align.csh align_batch.csh asa_cat asa_im_decode baseline_table.csh bperp calc_dop_orb \
	calc_dop_orb_envi cleanup.csh conv dem2topo_ra.csh detrend_before_unwrap.csh dump_orbit_envi.pl \
	dump_orbit_ers.pl dump_time_envi.pl ers_line_fixer esarp extend_orbit filter.csh find_auxi.pl \
	fitoffset.csh geocode.csh gmtsar.csh gmtsar_sharedir.csh grd2geotiff.csh grd2kml.csh intf.csh \
	intf_batch.csh landmask.csh make_a_offset.csh make_dem.csh make_los_ascii.csh make_profile.csh \
	make_raw_csk make_slc_csk make_slc_rs2 make_slc_s1a make_slc_tsx offset_topo p2p_ALOS.csh \
	p2p_ALOS2_SLC.csh p2p_ALOS_SLC.csh p2p_CSK.csh p2p_CSK_SLC.csh p2p_ENVI.csh p2p_ERS.csh \
	p2p_RS2_SLC.csh p2p_S1A_SLC.csh p2p_S1A_TOPS.csh p2p_SAT_SLC.csh p2p_TSX_SLC.csh phase2topo \
	phasediff phasefilt pre_proc.csh pre_proc_batch.csh pre_proc_init.csh proj_ll2ra.csh \
	proj_ll2ra_ascii.csh proj_model.csh proj_ra2ll.csh proj_ra2ll_ascii.csh read_data_file_ccrs \
	read_data_file_dpaf read_sarleader_dpaf resamp sarp.csh sbas slc2amp.csh snaphu snaphu.csh \
	snaphu_interp.csh stack.csh stack_corr.csh update_PRM.csh xcorr gmtsar_uninstall.sh"

printf "Remove: %s\n" $sharedir >&2
sudo rm -rf $sharedir
cd $dir
for exe in ${modules}; do
	printf "Remove: %s\n" $exe >&2
	sudo rm -f ${exe}
done
cd $here
