function [L1BS] = empty_L1BS_struct (L1BS)


global mode

L1BS = setfield(L1BS,'surf_counter',[]);
L1BS = setfield(L1BS,'lat_surf',[]);
L1BS = setfield(L1BS,'lon_surf',[]);
L1BS = setfield(L1BS,'alt_surf',[]);
L1BS = setfield(L1BS,'x_surf',[]);
L1BS = setfield(L1BS,'y_surf',[]);
L1BS = setfield(L1BS,'z_surf',[]);
L1BS = setfield(L1BS,'alt_rate_sat',[]);
L1BS = setfield(L1BS,'time_surf',[]);
L1BS = setfield(L1BS,'lat_sat',[]);
L1BS = setfield(L1BS,'lon_sat',[]);
L1BS = setfield(L1BS,'alt_sat',[]);
L1BS = setfield(L1BS,'x_sat',[]);
L1BS = setfield(L1BS,'y_sat',[]);
L1BS = setfield(L1BS,'z_sat',[]);
L1BS = setfield(L1BS,'x_vel_sat',[]);
L1BS = setfield(L1BS,'y_vel_sat',[]);
L1BS = setfield(L1BS,'z_vel_sat',[]);
L1BS = setfield(L1BS,'roll_surf',[]);
L1BS = setfield(L1BS,'pitch_surf',[]);
L1BS = setfield(L1BS,'yaw_surf',[]);
L1BS = setfield(L1BS,'win_delay_surf',[]);
L1BS = setfield(L1BS,'burst_index',[]);
L1BS = setfield(L1BS,'ProcessID_surf',[]);
L1BS = setfield(L1BS,'doppler_ang_surf',[]);
L1BS = setfield(L1BS,'beam_ang_surf',[]);
L1BS = setfield(L1BS,'ProcessID_beam',[]);
L1BS = setfield(L1BS,'x_vel_sat_beam',[]);
L1BS = setfield(L1BS,'y_vel_sat_beam',[]);
L1BS = setfield(L1BS,'z_vel_sat_beam',[]);
L1BS = setfield(L1BS,'x_sar_sat_beam',[]);
L1BS = setfield(L1BS,'y_sar_sat_beam',[]);
L1BS = setfield(L1BS,'z_sar_sat_beam',[]);
L1BS = setfield(L1BS,'win_delay_sar_ku_beam',[]);
L1BS = setfield(L1BS,'pri_sar_sat_beam',[]);

L1BS = setfield(L1BS,'Gap_flag',[]);
L1BS = setfield(L1BS,'pointing_ang_surf',[]);
L1BS = setfield(L1BS,'T0_sar_surf',[]);
L1BS = setfield(L1BS,'beams_surf',[]);
L1BS = setfield(L1BS,'N_beams_stack',[]);
L1BS = setfield(L1BS,'beam_index',[]);
L1BS = setfield(L1BS,'look_ang_surf',[]);

if(strcmp(mode,'SIN'))
    L1BS = setfield(L1BS,'beams_surf_2',[]);
end

% geometry corrections
L1BS = setfield(L1BS,'doppler_corr',[]);
L1BS = setfield(L1BS,'range_sat_surf',[]);
L1BS = setfield(L1BS,'slant_range_corr_time',[]);
L1BS = setfield(L1BS,'slant_range_corr',[]);
L1BS = setfield(L1BS,'wd_corr',[]);
L1BS = setfield(L1BS,'shift',[]);
L1BS = setfield(L1BS,'shift_coarse',[]);
L1BS = setfield(L1BS,'shift_fine',[]);
L1BS = setfield(L1BS,'N_windows',[]);
L1BS = setfield(L1BS,'good_samples',[]);
L1BS = setfield(L1BS,'beam_geo_corr',[]);

if(strcmp(mode,'SIN'))
    L1BS = setfield(L1BS,'beam_geo_corr_2',[]);
end
% range transformation

L1BS = setfield(L1BS,'beams_rng_cmpr',[]);
L1BS = setfield(L1BS,'beams_rng_cmprIQ',[]);

if(strcmp(mode,'SIN'))
    L1BS = setfield(L1BS,'beams_rng_cmpr_2' ,[]);
    L1BS = setfield(L1BS,'beams_rng_cmprIQ_2',[]);
    L1BS = setfield(L1BS,'phase_diff'       ,[]);
    L1BS = setfield(L1BS,'coherence'        ,[]);
    L1BS = setfield(L1BS,'coherence_mask'   ,[]);
end


% stack masking
L1BS = setfield(L1BS,'stack_mask',[]);


% multi-looking

% sigma0
L1BS = setfield(L1BS,'wfm_scaling_factor_sar_ku_beam',[]);



end

