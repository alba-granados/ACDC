% plrm_chain
function L1B = plrm_chain(L1A)
    global N_samples N_ku_pulses_burst_chd zp_fact_range_cnf
    global N_bursts_cycle_sar_chd pi_cst
    global chirp_slope_ku_chd wv_length_ku c_cst
    global power_tx_ant_ku_chd antenna_gain_ku_chd pulse_length_chd
    global earth_radius_cst sec_in_day_cst
    
    ref_burst = 2;
    wfm_L1A_fft_sq_av                   = zeros(N_bursts_cycle_sar_chd,N_samples*zp_fact_range_cnf);
    sigma0_scaling_factor_lrm_ku_burst  = zeros(1,N_bursts_cycle_sar_chd);
    altitude_rate                       = zeros(1,N_bursts_cycle_sar_chd);
    roll                                = zeros(1,N_bursts_cycle_sar_chd);
    pitch                               = zeros(1,N_bursts_cycle_sar_chd);
    yaw                                 = zeros(1,N_bursts_cycle_sar_chd);
    samples                             = 0:(N_samples-1);
    pulses_and_samples                  = repmat(samples,N_ku_pulses_burst_chd,1);
    
    for i_burst = 1:N_bursts_cycle_sar_chd
        %% Compute window delay shift
        wd_shift(i_burst) = (L1A(i_burst).win_delay_sar_ku - L1A(ref_burst).win_delay_sar_ku - (L1A(i_burst).alt_sar_sat - L1A(ref_burst).alt_sar_sat)*2/c_cst) ./ L1A(i_burst).T0_sar;
        
        %% Compute Doppler correction
        %(due to altitude rate) --> to add to the window delay
        doppler_seconds     = - 2 * L1A(i_burst).alt_rate_sar_sat / (wv_length_ku * chirp_slope_ku_chd); %in seconds
        doppler_meters      = doppler_seconds * c_cst / 2; %in meters
        
        %% Apply window delay shift
        wfm_L1A = squeeze(L1A(i_burst).wfm_cal_gain_corrected).* exp(2i*pi_cst/N_samples*wd_shift(i_burst).*pulses_and_samples);
        
        %% FFT range and power waveforms
        wfm_L1A_fft = abs(fftshift(fft(wfm_L1A.',N_samples*zp_fact_range_cnf),1).').^2/(N_samples);
%         wfm_L1A_fft = abs(fftshift(fft(wfm_L1A_fine.',N_samples*zp_fact_range_cnf),1).').^2/(N_samples*zp_fact_range_cnf);
        
        %% Intraburst AVERAGE
        wfm_L1A_fft_sq_av(i_burst,:) = mean(wfm_L1A_fft);
        
        %% Waveform scaling factor
        sigma0_scaling_factor_lrm_ku_burst(i_burst) =...
            10*log10(4) + 10*log10(pi_cst) - power_tx_ant_ku_chd...
            - 2*antenna_gain_ku_chd - 20*log10(wv_length_ku) + 10*log10(16)...
            + 10*log10(pi_cst) + 10*log10(pulse_length_chd * chirp_slope_ku_chd)...
            + 30*log10(L1A(i_burst).alt_sar_sat) - 10*log10(earth_radius_cst + L1A(i_burst).alt_sar_sat)...
            + 10*log10(earth_radius_cst) - 10*log10(c_cst);
        
        
        %% Parameters to average
        altitude_rate(i_burst)   = L1A(i_burst).alt_rate_sar_sat;
        roll(i_burst)            = L1A(i_burst).roll_sar;
        pitch(i_burst)           = L1A(i_burst).pitch_sar;
    	yaw(i_burst)             = L1A(i_burst).yaw_sar;
        
    end
    
%     toc(t0)


    %% Window delay and altitude adjustment --> all bursts now have the same value
    L1B.win_delay_sar           = L1A(ref_burst).win_delay_sar_ku;
    L1B.altitude                = L1A(ref_burst).alt_sar_sat;
    
    %% Interburst AVERAGE
    L1B.wfm_L1A_fft_sq_av_RC        = mean(wfm_L1A_fft_sq_av);
    L1B.sigma0_scaling_factor_RC    = mean(sigma0_scaling_factor_lrm_ku_burst);
    
    
    
    
    %% Prepare Writing
    %---------- A. Time variables -----------------------------------------
    L1B.time_l1b_plrm               = L1A(ref_burst).days + L1A(ref_burst).seconds + L1A(ref_burst).microseconds; %to check units of microseconds
    
    %---------- B. Orbit and attitude variables ---------------------------
    L1B.latitude                    = L1A(ref_burst).lat_sar_sat;
    L1B.longitude                   = L1A(ref_burst).lon_sar_sat;
    if (L1B.longitude > 180)
        L1B.longitude               = L1B.longitude - 360; 
    end
    
    L1B.altitude_rate               = mean(altitude_rate);
    
    L1B.roll                        = mean(roll);
    L1B.pitch                       = mean(pitch);
    L1B.yaw                         = mean(yaw);
    
    %----------D. Altimeter range variables -------------------------------
    L1B.doppler                     = mean(doppler_meters);
    L1B.win_delay_sar               = L1B.win_delay_sar + L1B.doppler;
    
    %----------L. Altimeter engineering variables -------------------------
    L1B.T0_sar                      = L1A(ref_burst).T0_sar;
    
    
% end







% % % %%
% % % % if WRITE
% % %     %*****************************************************************************************
% % %     %****************************** L1B WRITING **********************************************
% % %     %*****************************************************************************************
% % % %     if strcmp(TEST,'T01')
% % % %         filename = 'S6_P4_SIM_LRM_L1B__20190119T064000_20190119T064019_T01.nc';
% % % %     elseif strcmp(TEST,'T02')
% % % %         filename = 'S6_P4_SIM_LRM_L1B__20190119T064000_20190119T064019_T01.nc';
% % % %     elseif strcmp(TEST,'T03')
% % % %         filename = 'S6_P4_SIM_LRM_L1B__20190119T064000_20190119T064019_T01.nc';
% % % %     end
% % %     
% % %     
% % %     %---------- Preparation ---------------------------------------------------
% % %     dimensions_key = 'Dimensions'; format_key = 'Format'; data_type_key = 'DataType';
% % %     
% % %     netcdf_v4_format = 'netcdf4';
% % %     
% % %     ku_rec_dimension = 'Ku_rec'; np_dimension = 'Np'; ns_dimension = 'Ns';
% % %     space_3D_dimension = 'space_3D'; space_3D_dimension_size = 3;
% % %     
% % %     long_name_att = 'long_name'; comment_att = 'comment';
% % %     units_att = 'units'; scale_factor_att = 'scale_factor';
% % %     add_offset_att = 'add_offset';
% % %     
% % %     int8_type = 'int8'; uint8_type = 'uint8'; int16_type = 'int16';
% % %     uint16_type = 'uint16'; int32_type = 'int32'; uint32_type = 'uint32';
% % %     uint64_type = 'uint64'; float_type = 'single';
% % %     
% % %     day_units = 'day'; seconds_units = 'seconds';
% % %     number_units = '#'; degrees_units = 'degrees';
% % %     meters_units = 'meters'; meters_per_second_units = 'm/s';
% % %     dB_units = 'dB'; Hz_units = 'Hz'; T0x4_units = 'T0*4';
% % %     T0d64_units = 'T0/64'; T0x8_units = 'T0*8'; T0d16d64_units = 'T0/16/64';
% % %     sqrtW_per_count_units = 'sqrt(Watt)/#';
% % %     
% % %     
% % %     
% % %     %---------- A. Time variables ---------------------------------------------
% % %     l1_mode_id_ku = L1A.data.l1_mode_id_ku(1:N_bursts_cycle_sar_chd:N_bursts_cycle_sar_chd);
% % %     % time_day_lrm - Done above
% % %     % time_seconds_lrm - Done above
% % %     tm_source_sequence_counter_ku = L1A.data.tm_source_sequence_counter_ku(1:N_bursts_cycle_sar_chd:N_bursts_cycle_sar_chd);
% % %     l1b_record_counter_ku = L1A.data.l1b_record_counter_ku(1:N_bursts_cycle_sar_chd:N_bursts_cycle_sar_chd);
% % %     
% % %     
% % %     %---------- B. Orbit and attitude variables -------------------------------
% % %     latitude_ku = int32(lat_sat * 1e7);
% % %     for i_lon = 1:size(lon_sat,1)
% % %         if (lon_sat(i_lon) > 180)
% % %             longitude_ku = uint32((lon_sat +360)* 1e7); 
% % %         else
% % %             longitude_ku = uint32((lon_sat )* 1e7); 
% % %         end
% % %     end
% % %     com_altitude_ku = uint32((alt_sat - 1.3e6) * 1e4);
% % %     com_altitude_rate_ku = int32(alt_rate_sat * 1e4);
% % %     com_velocity_vector_ku = int32([x_vel_sat * 1e4; y_vel_sat * 1e4; z_vel_sat * 1e4]);
% % %     satellite_mispointing_ku = int32([pitch_surf * 180/pi_cst * 1e7; roll_surf * 180/pi_cst * 1e7; yaw_surf * 180/pi_cst * 1e7]);
% % %     pitch_bias = 0; % to CHD file
% % %     roll_bias = 0; % to CHD file
% % %     yaw_bias = 0; % to CHD file
% % %     mispointing_bias_ku = int32([pitch_bias * 180/pi_cst * 1e7, roll_bias * 180/pi_cst * 1e7, yaw_bias * 180/pi_cst * 1e7]);
% % %     
% % %     
% % %     %----------C. Configuration and quality variables -------------------------
% % %     l1_instrument_configuration_ku = zeros(N_total_bursts_sar_ku_isp,1);
% % %     for i_burst=1:N_total_bursts_sar_ku_isp
% % %         bit0    = dec2bin(inst_id_sar_isp(i_burst));
% % %         bit15   = dec2bin(tm_mode_id_sar_isp(i_burst),4);
% % %         trk_bin = dec2bin(trk_config_sar_isp(i_burst),8);
% % %         bit68   = trk_bin((1:3));
% % %         bit9    = dec2bin(loss_track_criterion_isp(i_burst));
% % %         bit10   = dec2bin(nav_bulletin_status_isp(i_burst));
% % %         bit11to15= '0000';
% % %          l1_instrument_configuration_ku(i_burst,1) = bin2dec([bit10 bit9 bit68 bit15 bit0]);
% % %     %     l1_instrument_configuration_ku(i_burst,1) = bin2dec([bit11to15 bit10 bit9 bit68 bit15 bit0]);
% % %     %    l1_instrument_configuration_ku(i_burst,1) = bin2dec([bit0 bit15 bit68 bit9 bit10 bit11to15]);
% % %     end
% % %     
% % %     for i_surf = 1:N_total_surf_loc
% % %         l1_instrument_configuration_ku_surf(i_surf)     = l1_instrument_configuration_ku(burst_index_nadir(i_surf));
% % %     end
% % %     
% % %     l1b_mcd_ku = zeros(N_total_surf_loc,1); % TBD
% % %     
% % %     
% % %     %----------D. Altimeter range variables -----------------------------------
% % %     altimeter_range_calibrated_ku = uint32(((win_delay_surf) * c_cst/2 - 1.3e6) * 1e4);
% % %     % altimeter_range_calibrated_ku = uint32(((zeros(1,N_total_surf_loc)+win_delay_surf(1)) * c_cst/2 - 1.3e6) * 1e4);
% % %     range_corr_internal_delay_ku = int16(int_delay_cor_cal1_surf * c_cst/2 * 1e4);
% % %     range_corr_com_ku = int16(-(z_cog_ant_corr_surf.') * 1e4);
% % %     
% % %     
% % %     %----------E. Altimeter power variables -----------------------------------
% % %     attenuator_calibrated_ku = int16(att_conv_sar_ku_surf * 1e2);
% % %     altimeter_power_drift_ku  = int16(power_var_cal1_sar_ku_surf * 1e2);
% % %     power_corr_digital_processing_ku = int16(onboard_proc_sar * 1e2);
% % %     power_scaling_to_antenna_ku = int16(gain_corr_instr_sar_surf * 1e2);
% % %     
% % %     
% % %     %----------F. Altimeter engineering variables -----------------------------
% % %     altimeter_clock_ku = int32(1./T0_sar_surf_nadir - 3.95e8)* 1e9;
% % %     tm_h0_ku = uint32(h0_sar_isp_surf);
% % %     tm_cor2_ku = int16(cor2_sar_isp_surf);
% % %     hn_mean_ku = uint32(mean_cai_fai_sar_surf);
% % %     pri_lrm_l1b_ku = uint32(pri_sar_isp_surf .* T0_sar_surf_nadir * pri_T0_unit_conv_chd * 1e12);
% % %     
% % %     
% % %     %----------M. Waveform related variables ----------------------------------
% % %     waveform_scale_factor_ku = zeros(N_total_surf_loc,1);
% % %     sar_power_waveform_ku = zeros(N_total_surf_loc, N_samples_sar_chd * zp_fact_range_cnf);
% % %     
% % %     if(no_amb==0)
% % %         for i_surf=1:N_total_surf_loc
% % %             waveform_scale_factor_ku(i_surf) = single(max(wfm_cor_i2q2_sar_ku(i_surf,:))*10^0.3 / (2^16-1));
% % %             sar_power_waveform_ku(i_surf,:) = uint16(round(wfm_cor_i2q2_sar_ku(i_surf,:).*10^0.3 ./ waveform_scale_factor_ku(i_surf)));
% % %         end
% % %     elseif(no_amb==1)
% % %        for i_surf=1:N_total_surf_loc
% % %             waveform_scale_factor_ku(i_surf) = single(max(wfm_cor_i2q2_sar_ku_noDopp(i_surf,:))*10^0.3 / (2^16-1));
% % %             sar_power_waveform_ku(i_surf,:) = uint16(round(wfm_cor_i2q2_sar_ku_noDopp(i_surf,:).*10^0.3 ./ waveform_scale_factor_ku(i_surf)));
% % %         end     
% % %     end
% % %     snr_estimation_ku = zeros(N_total_surf_loc,1); % TBD
% % %     sigma0_scaling_factor_ku = int16(wfm_scaling_factor_sar_ku * 1e2);
% % %     
% % %     
% % %     
% % % % end








