function [burst] = read_FBR_CryoSat_time_orbit_group(fid, burst)

global sec_in_day_cst pri_sar_chd N_bursts_cycle_chd

    %--------------------------------%
    %--  read Time and Orbit Group --%
    %--------------------------------%
      
        
		%1 Data Record Time (MDSR TimeStamp) TAI 12 sl+2*ul
        burst.days = fread(fid,1,'int32')* sec_in_day_cst;
        
        if(isempty(burst.days)==1)%If we don't find anything inside means that the file has end
            end_of_file=1;
            record_num = record_num-1;
        else
        
            burst.seconds = fread(fid,1,'uint32');
            burst.microseconds = fread(fid, 1, 'uint32')* 1e-6;
            
            %             data_record_time = struct(  'day', day,...
%                                         'second', second, ...
%                                         'microsecond',microsecond); 
           
            %2 USO Correction Factor (ratio) 10-9 4 sl
            burst.USO_correction = fread(fid,1,'int32')*1e-15;
            %3 Mode ID 2 us (see table 2.3.3-2) 
            modeID_mode = fread(fid,1,'ubit6');
            modeID_SARIn_deg_case = fread(fid,1,'ubit1');
            modeID_reserved = fread(fid,1,'ubit1');
            modeID_CAL4 = fread(fid,1,'ubit1');
            modeID_Platform_Cont = fread(fid,1,'ubit2');
            modeID_reserved2 = fread(fid,1,'ubit5'); 

%             modeID_CR = struct('modeID_mode', modeID_mode,...
%                                     'modeID_SARIn_deg_case', modeID_SARIn_deg_case,...
%                                     'modeID_reserved', modeID_reserved,...
%                                     'modeID_CAL4', modeID_CAL4,...
%                                     'modeID_Platform_Cont', modeID_Platform_Cont,...
%                                     'modeID_reserved2',modeID_reserved2);
           switch modeID_mode
              case 2 %SAR
                   burst.ProcessID = 58;
              case 3 %SARIN
                   burst.ProcessID = 58;
                   if modeID_CAL4
                       
                       burst.ProcessID= 57; %CAL4 burst
                       
                   end
               otherwise
                   burst.ProcessID = 0;
                       
                   
           end
                   
            
            %4 Source Sequence Counter 2 us (see note 6)
            burst.source_seq_count_sar_ku_fbr = fread(fid,1,'uint16');
            %burst.source_seq_count_sar_ku_fbr
            %5 Instrument Configuration 4 ul (see table 2.3.3-3)

            ins_rx_in_use = fread(fid,1,'ubit2');
            burst.ins_id = fread(fid,1,'ubit1');
            ins_reserved = fread(fid,1,'ubit1');
            ins_bandwidth = fread(fid,1,'ubit2');
            ins_reserved2 = fread(fid,1,'ubit2');
            burst.ins_tracking_mode = fread(fid,1,'ubit2');
            ins_ext_calibration = fread(fid,1,'ubit1');
            ins_reserved3 = fread(fid,1,'ubit1');
            burst.ins_loop_stat = fread(fid,1,'ubit1');
            ins_loss_echo = fread(fid,1,'ubit1');
            ins_real_time_error = fread(fid,1,'ubit1');
            ins_echo_sat_error = fread(fid,1,'ubit1');
            ins_rx_band_attenuation = fread(fid,1,'ubit1');
            ins_cycle_report = fread(fid,1,'ubit1');
            ins_star_tracker_1 = fread(fid,1,'ubit1');
            ins_star_tracker_2 = fread(fid,1,'ubit1');
            ins_star_tracker_3 = fread(fid,1,'ubit1');
            ins_reserved4 = fread(fid,1,'ubit11');
            
            burst.inst_id_sar_isp=0;


            %6 Burst counter (always starts from 1 and incremented at group rate) 4 ul
            
            burst.pri_sar_isp = pri_sar_chd;
            burst.ambiguity_order_sar_isp = 0;
            burst.burst_sar_ku = fread(fid,1,'uint32');
            %burst.burst_sar_ku
            burst.burst_sar_ku_fbr= 1+mod(burst.burst_sar_ku,N_bursts_cycle_chd);
%             if(burst.burst_sar_ku_fbr==0)
%                 burst.burst_sar_ku_fbr=1;
%             end
            %7 Latitude of measurement 10-1 mcrodeg 4 sl (see note 1)
            burst.lat_sar_sat = fread (fid,1,'int32') * 1e-7;
            
            %8 Longitude of measurement 10-1 microdeg 4 sl (see note 1)
            burst.lon_sar_sat = fread (fid,1,'int32') * 1e-7;
            
            %9 Altitude of COG above reference ellipsoid (interpolated value)mm 4 sl
            burst.alt_sar_sat = fread (fid,1,'int32') * 1e-3;
            
            %10 Instantaneous altitude rate derived from orbit mm/s 4 sl
            burst.alt_rate_sar_sat = fread (fid,1,'int32') * 1e-3;
            
            
            %11 Satellite velocity vector[3](in ITRF) mm/s 3*4 sl
            burst.x_vel_sat_sar = fread(fid,1,'int32') * 1e-3;
            
            burst.y_vel_sat_sar = fread(fid,1,'int32') * 1e-3;
            
            burst.z_vel_sat_sar = fread(fid,1,'int32') * 1e-3;
            
            
            %12 Real beam direction vector[3](in CRF) micros 3*4 sl
            real_beam_direction(1) = fread(fid,1,'int32')* 1e-6;
            real_beam_direction(2) = fread(fid,1,'int32')* 1e-6;
            real_beam_direction(3) = fread(fid,1,'int32')* 1e-6;
            %13 Interferometer baseline vector[3](in CRF)micros 3*4 sl
            inferometer_baseline(1) = fread(fid,1,'int32')* 1e-6;
            inferometer_baseline(2) = fread(fid,1,'int32')* 1e-6;
            inferometer_baseline(3) = fread(fid,1,'int32')* 1e-6;
            
            burst.roll_sar  = atan((-1.0*inferometer_baseline(1)/inferometer_baseline(3))); % [rad] inferometer_baseline(1);%
            burst.pitch_sar = atan((-real_beam_direction(2)/real_beam_direction(1))); % [rad]-real_beam_direction(2);
            burst.yaw_sar   = atan((inferometer_baseline(2)/inferometer_baseline(3))); % [rad]-inferometer_baseline(2);%
            
        
            
            %14 FBR Measurement ConfidenceData (flag word)4 u
            burst.confi_block_degraded = fread(fid,1,'ubit1');
                confi_blank_block = fread(fid,1,'ubit1');
                confi_datation_degraded = fread(fid,1,'ubit1');
                confi_orbit_prop_error = fread(fid,1,'ubit1');
                confi_orbit_file_change = fread(fid,1,'ubit1');
                confi_orbit_discon = fread(fid,1,'ubit1');
                confi_echo_sat = fread(fid,1,'ubit1');
                confi_other_echo_error = fread(fid,1,'ubit1');
                confi_rx1_error_for_SARIN = fread(fid,1,'ubit1');
                confi_rx2_error_for_SARIN = fread(fid,1,'ubit1');
                confi_window_delay_inconsistency = fread(fid,1,'ubit1');
                confi_AGC_inconsistency = fread(fid,1,'ubit1');
                confi_cal1_correction_miss = fread(fid,1,'ubit1');
                confi_cal1_correction_from_IPF_DB = fread(fid,1,'ubit1');
                confi_DORIS_USO_correction = fread(fid,1,'ubit1');
                confi_complex_cal1_correction_from_IPF_DB = fread(fid,1,'ubit1');
                confi_TRK_echo_error = fread(fid,1,'ubit1');
                confi_echo_rx1_error = fread(fid,1,'ubit1');
                confi_echo_rx2_error = fread(fid,1,'ubit1');
                confi_NMP_inconsistency = fread(fid,1,'ubit1');
                confi_azimuth_cal_missing = fread(fid,1,'ubit1');
                confi_azimuth_cal_from_IPF_DB = fread(fid,1,'ubit1');
                confi_range_window_cal_function_missing = fread(fid,1,'ubit1');
                confi_range_window_cal_function_from_IPF_DB = fread(fid,1,'ubit1');
                confi_reserved = fread(fid,1,'ubit1');
                confi_cal2_correction_missing = fread(fid,1,'ubit1');
                confi_cal2_correction_from_IPF_DB = fread(fid,1,'ubit1');
                confi_power_scaling_error_LRM_only = fread(fid,1,'ubit1');
                confi_attitude_correction_missing = fread(fid,1,'ubit1');
                confi_attitude_interpolation_error = fread(fid,1,'ubit1');
                confi_reserved2 = fread(fid,1,'ubit1');
                confi_phase_perturbation = fread(fid,1,'ubit1');
    
%                 mea_conf_data_sar_ku_fbr=0; %TBD
%                 burst.mea_conf_data_sar_ku_fbr =struct(    'confi_block_degraded', burst.confi_block_degraded,...
%                                     'confi_blank_block', confi_blank_block,...
%                                     'confi_datation_degraded', confi_datation_degraded,...
%                                     'confi_orbit_prop_error', confi_orbit_prop_error,...
%                                     'confi_orbit_file_change', confi_orbit_file_change,...
%                                     'confi_orbit_discon',confi_orbit_discon,...                                
%                                     'confi_echo_sat', confi_echo_sat,...
%                                     'confi_other_echo_error', confi_other_echo_error,...
%                                     'confi_rx1_error_for_SARIN', confi_rx1_error_for_SARIN,...
%                                     'confi_rx2_error_for_SARIN', confi_rx2_error_for_SARIN,...
%                                     'confi_window_delay_inconsistency',confi_window_delay_inconsistency,...
%                                     'confi_AGC_inconsistency', confi_AGC_inconsistency,...
%                                     'confi_cal1_correction_miss', confi_cal1_correction_miss,...
%                                     'confi_cal1_correction_from_IPF_DB', confi_cal1_correction_from_IPF_DB,...
%                                     'confi_DORIS_USO_correction', confi_DORIS_USO_correction,...
%                                     'confi_complex_cal1_correction_from_IPF_DB', ins_star_tracker_2,...
%                                     'confi_TRK_echo_error', ins_star_tracker_3,...
%                                     'confi_echo_rx1_error',ins_reserved4,...                            
%                                     'confi_echo_rx2_error', confi_echo_rx2_error,...
%                                     'confi_NMP_inconsistency', confi_NMP_inconsistency,...
%                                     'confi_azimuth_cal_missing', confi_azimuth_cal_missing,...
%                                     'confi_azimuth_cal_from_IPF_DB', confi_azimuth_cal_from_IPF_DB,...
%                                     'confi_range_window_cal_function_missing',confi_range_window_cal_function_missing,...
%                                     'confi_range_window_cal_function_from_IPF_DB', confi_range_window_cal_function_from_IPF_DB,...
%                                     'confi_reserved', confi_reserved,...
%                                     'confi_cal2_correction_missing', confi_cal2_correction_missing,...
%                                     'confi_cal2_correction_from_IPF_DB', confi_cal2_correction_from_IPF_DB,...
%                                     'confi_power_scaling_error_LRM_only', confi_power_scaling_error_LRM_only,...
%                                     'confi_attitude_correction_missing', confi_attitude_correction_missing,...
%                                     'confi_attitude_interpolation_error',confi_attitude_interpolation_error,...
%                                     'confi_reserved2',confi_reserved2,...
%                                     'confi_phase_perturbation',confi_phase_perturbation);


    
        end
end