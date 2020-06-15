function [L1A] = empty_L1A_struct(L1A)


% global N_ku_pulses_burst_chd N_samples_sar_chd
global mode

for s=1:length(L1A)
    L1A(s) = setfield_new(L1A(s),'source_seq_count_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'days',[]);
    L1A(s) = setfield_new(L1A(s),'seconds',[]);
    L1A(s) = setfield_new(L1A(s),'microseconds',[]);
    L1A(s) = setfield_new(L1A(s),'USO_correction',[]);
    L1A(s) = setfield_new(L1A(s),'ProcessID',[]);
    L1A(s) = setfield_new(L1A(s),'source_seq_count_sar_ku_fbr',[]);
    L1A(s) = setfield_new(L1A(s),'ins_id',[]);
    L1A(s) = setfield_new(L1A(s),'ins_tracking_mode',[]);
    L1A(s) = setfield_new(L1A(s),'ins_loop_stat',[]);
    L1A(s) = setfield_new(L1A(s),'inst_id_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'pri_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'ambiguity_order_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'burst_sar_ku',[]);
    L1A(s) = setfield_new(L1A(s),'burst_sar_ku_fbr',[]);
    L1A(s) = setfield_new(L1A(s),'lat_sar_sat',[]);
    L1A(s) = setfield_new(L1A(s),'lon_sar_sat',[]);
    L1A(s) = setfield_new(L1A(s),'alt_sar_sat',[]);
    L1A(s) = setfield_new(L1A(s),'alt_rate_sar_sat',[]);
    L1A(s) = setfield_new(L1A(s),'x_vel_sat_sar',[]);
    L1A(s) = setfield_new(L1A(s),'y_vel_sat_sar',[]);
    L1A(s) = setfield_new(L1A(s),'z_vel_sat_sar',[]);
    L1A(s) = setfield_new(L1A(s),'roll_sar',[]);
    L1A(s) = setfield_new(L1A(s),'pitch_sar',[]);
    L1A(s) = setfield_new(L1A(s),'yaw_sar',[]);
    L1A(s) = setfield_new(L1A(s),'confi_block_degraded',[]);
    L1A(s) = setfield_new(L1A(s),'mea_conf_data_sar_ku_fbr',[]);
    L1A(s) = setfield_new(L1A(s),'win_delay_sar_ku',[]);
    if(strcmp(mode,'SIN'))
        L1A(s) = setfield_new(L1A(s),'win_delay_sar_ku_2',[]);
    end
    L1A(s) = setfield_new(L1A(s),'h0_comp_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'cor2_comp_sar_isp',[]);
    %L1A(s) = setfield_new(L1A(s),'h0_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'cai_sar_isp',[]);
    %L1A(s) = setfield_new(L1A(s),'cor2_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'fai_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'ATT1_science',[]);
    L1A(s) = setfield_new(L1A(s),'ATT2_science',[]);
    L1A(s) = setfield_new(L1A(s),'att_sar_ku_isp',[]);
    L1A(s) = setfield_new(L1A(s),'tot_fixed_gain_1',[]);
    L1A(s) = setfield_new(L1A(s),'tot_fixed_gain_2',[]);
    L1A(s) = setfield_new(L1A(s),'transmit_power',[]);
    L1A(s) = setfield_new(L1A(s),'doppler_range_correction',[]);
    L1A(s) = setfield_new(L1A(s),'instrument_range_correction_tx_rx',[]);
    L1A(s) = setfield_new(L1A(s),'instrument_range_correction_rx',[]);
    L1A(s) = setfield_new(L1A(s),'instrument_sigma0_correction_tx_rx',[]);
    L1A(s) = setfield_new(L1A(s),'instrument_sigma0_correction_rx',[]);
    L1A(s) = setfield_new(L1A(s),'internal_phase_correction',[]);
    L1A(s) = setfield_new(L1A(s),'external_phase_correction',[]);
    L1A(s) = setfield_new(L1A(s),'noise_power',[]);
    L1A(s) = setfield_new(L1A(s),'phase_slope_correction',[]);
    L1A(s) = setfield_new(L1A(s),'gain_corr_instr_sar_1',[]);
    if(strcmp(mode,'SIN'))
        L1A(s) = setfield_new(L1A(s),'gain_corr_instr_sar_2',[]);
    end


    L1A(s) = setfield_new(L1A(s),'dry_tropo_correction_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'wet_tropo_correction_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'inverse_baro_correction_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'Dynamic_atmospheric_correction_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'GIM_iono_correction_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'model_iono_correction_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'ocean_equilibrium_tide_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'long_period_tide_height_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'ocean_loading_tide_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'solid_earth_tide_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'geocentric_polar_tide_bursts',[]);
    L1A(s) = setfield_new(L1A(s),'surface_type_flag_bursts',[]);

    L1A(s) = setfield_new(L1A(s),'nimp_sar_isp',[]);
    L1A(s) = setfield_new(L1A(s),'wfm_cal_gain_corrected',[]);
    if(strcmp(mode,'SIN'))
        L1A(s) = setfield_new(L1A(s),'wfm_cal_gain_corrected_2',[]);
    end

    L1A(s) = setfield_new(L1A(s),'time_sar_ku',[]);
    L1A(s) = setfield_new(L1A(s),'x_sar_sat',[]);
    L1A(s) = setfield_new(L1A(s),'y_sar_sat',[]);
    L1A(s) = setfield_new(L1A(s),'z_sar_sat',[]);
    L1A(s) = setfield_new(L1A(s),'lat_sar_surf',[]);
    L1A(s) = setfield_new(L1A(s),'lon_sar_surf',[]);
    L1A(s) = setfield_new(L1A(s),'alt_sar_surf',[]);
    L1A(s) = setfield_new(L1A(s),'x_sar_surf',[]);
    L1A(s) = setfield_new(L1A(s),'y_sar_surf',[]);
    L1A(s) = setfield_new(L1A(s),'z_sar_surf',[]);
    L1A(s) = setfield_new(L1A(s),'doppler_ang_sar_sat',[]);

    L1A(s) = setfield_new(L1A(s),'T0_sar',[]);
    L1A(s) = setfield_new(L1A(s),'pri_sar',[]);

    %beam angles
    L1A(s) = setfield_new(L1A(s),'N_beams_sar_ku',0);
    L1A(s) = setfield_new(L1A(s),'surf_loc_index',0);
    L1A(s) = setfield_new(L1A(s),'beam_ang',[]);
    L1A(s) = setfield_new(L1A(s),'beam_ang_nadir_index',[]);
    L1A(s) = setfield_new(L1A(s),'beam_ang_index',[]);
    L1A(s) = setfield_new(L1A(s),'start_beam',[]);
    L1A(s) = setfield_new(L1A(s),'end_beam',[]);
    L1A(s) = setfield_new(L1A(s),'beam_index',[]);

    %azimuth processing
    L1A(s) = setfield_new(L1A(s),'beams_focused_shifted',[]);
    if(strcmp(mode,'SIN'))
        L1A(s) = setfield_new(L1A(s),'beams_focused_shifted_2',[]);
    end

%     % added for FAI validation purposes
%     L1A(s) = setfield_new(L1A(s),'wfm_cal_gain_corrected_before_align',[]);
%     L1A(s) = setfield_new(L1A(s),'wfm_cal_gain_corrected_FAI_positive_sign',[]);    
%     L1A(s) = setfield_new(L1A(s),'idx_leading_noFAI',[]);
%     L1A(s) = setfield_new(L1A(s),'idx_leading_FAI',[]);
%     L1A(s) = setfield_new(L1A(s),'idx_leading_FAI_plussign',[]);
    
end
end

