configure_file (gmtsar_sharedir.csh.in gmtsar_sharedir.csh)

install (PROGRAMS align.csh align_ALOS2_SCAN.csh align_ALOS_SLC.csh align_batch.csh align_batch_ALOS2_SCAN.csh
	align_batch_ALOS_SLC.csh align_tops.csh align_tops_6par.csh align_tops_esd.csh baseline_table.csh
	cleanup.csh create_frame_tops.csh dem2topo_ra.csh dem2topo_ra_ALOS2.csh detrend_before_unwrap.csh
	filter.csh fitoffset.csh geocode.csh gmtsar.csh gmtsar_sharedir.csh grd2geotiff.csh grd2kml.csh intf.csh
	intf_batch.csh intf_batch_ALOS2_SCAN.csh intf_tops.csh landmask.csh landmask_ALOS2.csh m2s.csh
	make_a_offset.csh make_dem.csh make_los_ascii.csh make_profile.csh merge_batch.csh
    merge_unwrap_geocode_tops.csh p2p_ALOS2_SCAN_Frame.csh p2p_ALOS2_SCAN_SLC.csh
	p2p_ENVI.csh p2p_ERS.csh p2p_S1_TOPS.csh p2p_S1_TOPS_Frame.csh pre_proc.csh pre_proc_batch.csh
	pre_proc_batch_ALOS2_SCAN.csh pre_proc_batch_ALOS_SLC.csh pre_proc_init.csh preproc_batch_tops.csh
	preproc_batch_tops_esd.csh proj_ll2ra.csh proj_ll2ra_ascii.csh proj_model.csh proj_ra2ll.csh
	proj_ra2ll_ascii.csh sarp.csh slc2amp.csh snaphu.csh snaphu_interp.csh stack.csh stack_corr.csh 
	update_PRM.csh select_pairs.csh
	DESTINATION bin
	COMPONENT Runtime)
