% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
% 
% This code implements the algorithm as described in the
% ISARD_ESA_CR2_TRP_CAL_DPM_030 2.b of 26/05/2011
%
% ---------------------------------------------------------
% READFBR: function that reads the FBR data set from input filename,
% assuming SARIN data records
%
% Calling
%   fbr_ds = readFBR( filename, headers )
%
% Inputs
%   filename: input SARIN FBR file
%   headers:   if false --> without header
%              if true --> with header
%
% Output
%   fbr_ds           : data contained in the file 
%
% ----------------------------------------------------------
% 
% Author:   Albert Garcia / isardSAT
%           Eduard Makhoul / isardSAT
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Mònica Roca / isardSAT (26/05/11)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%modify the inputs to provide optionally the filenames of
%cal1-pulse-to-pulse and CAL2-
function [data] = readFBR_CryoSat_SAR(filename_FBR,first_record, num_records)


global N_ku_pulses_burst_chd N_samples_sar_chd sec_in_day_cst N_bursts_cycle_chd
global flat_coeff_cst semi_major_axis_cst 
global pri_sar_chd bri_chd prf_sar_chd brf_chd T0_chd
global pri_sar_nom bri_nom prf_sar_nom brf_nom c_cst

%modified by EM 30.03.2016: use global variables to define the
%CAl1-pulse-to-pulse and CAL2 corrections
global burst_power_array_cor_cal1_sar_rep burst_phase_array_cor_cal1_sar_rep wfm_cal2_science_sar_rep 

t1 = tic;
% filename = 'input_data/test converter/CS_OFFL_SIR1SAR_FR_20110910T140249_20110910T140723_B001.DBL';
end_of_file = 0;
shift_isp = 0;
record_num=0;


    disp(['reading FBR from ' filename_FBR ' WITH headers']);
    sph = readSPH(filename_FBR);
    % FBR data set must be the first one
    fbr_dsd = sph.dsds(1);
    if deblank(fbr_dsd.ds_name) ~= 'SIR_FBR_SAR'
        disp('no SIR_FBR_SAR data set found');
        return;
    else
        disp('found SIR_FBR_SARIN data set');
        disp(fbr_dsd);
    end

    ds_offset = fbr_dsd.ds_offset;
    ds_size = fbr_dsd.ds_size;
%      num_dsr = fbr_dsd.num_dsr;
    dsr_size = fbr_dsd.dsr_size;

     num_dsr=num_records;
     

% read the dsrs
fid = fopen(filename_FBR,'r','b');
fseek(fid,ds_offset+(first_record-1)*dsr_size,'bof'); % last records of the product
% fseek(fid,ds_offset,'bof'); % start reading from the first record

end_of_file = 0;
i_burst = 0;

% modes (by modeID_record, 32bits, TABLE02)
LRM_KU = 0;
LRM_C = 1024;
SAR_KU = 2048;
SAR_C = 3072;
CAL1_LRM_I2Q2_KU = 10240;
CAL1_LRM_I2Q2_C = 11264;
CAL1_LRM_IQ_KU = 12288;
CAL1_LRM_IQ_C = 13312;
CAL1_SAR_KU = 14336;
CAL1_SAR_C = 15360;
 
 wfm_iq_sar_ku_fbr_aux_1 = zeros(1,N_ku_pulses_burst_chd * N_samples_sar_chd*2);
%  wfm_iq_sar_ku_fbr_aux_2 = zeros(1,N_ku_pulses_burst_chd * N_samples_sar_chd*2);
 
 wfm_iq_sar_ku_fbr_i_11     = zeros(num_dsr*20,N_ku_pulses_burst_chd ,N_samples_sar_chd);
 wfm_iq_sar_ku_fbr_q_11     = zeros(num_dsr*20,N_ku_pulses_burst_chd ,N_samples_sar_chd);

%  wfm_iq_sar_ku_fbr_i_22     = zeros(num_dsr*20,N_ku_pulses_burst_chd ,N_samples_sar_chd);
%  wfm_iq_sar_ku_fbr_q_22     = zeros(num_dsr*20,N_ku_pulses_burst_chd ,N_samples_sar_chd);

%  wfm_iq_sar_ku_fbr_i_1 = zeros(num_dsr*20,N_ku_pulses_burst_chd,N_samples_sar_chd);
%  wfm_iq_sar_ku_fbr_q_1 = zeros(num_dsr*20,N_ku_pulses_burst_chd,N_samples_sar_chd); 
%  wfm_iq_sar_ku_fbr_i_2 = zeros(num_dsr*20,N_ku_pulses_burst_chd,N_samples_sar_chd);
%  wfm_iq_sar_ku_fbr_q_2 = zeros(num_dsr*20,N_ku_pulses_burst_chd,N_samples_sar_chd);


% progressbar('Records','Data Blocks')
%Added by EM 30.03.2016
% Flags to activate or not CAL1-pulse-to-pulse and CAL2 based on the values
% in the corresponding variables 
%CAL1-pulse-to-pulse
CAL1pp_active=~isempty(find(burst_power_array_cor_cal1_sar_rep, 1)) || ~isempty(find(burst_phase_array_cor_cal1_sar_rep, 1));
%CAL2
CAL2_active=~isempty(find(~wfm_cal2_science_sar_rep, 1));

%% READ FBR FILE 
for record_num=1:num_dsr
    
    
   
    
    %--------------------------------%
    %--  read Time and Orbit Group --%
    %--------------------------------%
    j = 1;
    
    
    for j=1:20
%         progressbar('Read FBR','Surface Locations',...
%         'Beam Angles','Azimuth Processing','Stacking','Geometry Corrections',...
%         'ACDC','Range Compression','Multi-Looking');

        i_burst= (record_num-1)*20+j;
         
        data.source_seq_count_sar_isp(i_burst)=record_num;
        
		%1 Data Record Time (MDSR TimeStamp) TAI 12 sl+2*ul
        data.days(i_burst) = fread(fid,1,'int32')* sec_in_day_cst;
        
        if(isempty(data.days(i_burst))==1)%If we don't find anything inside means that the file has end
            end_of_file=1;
            record_num = record_num-1;
        else
        
            data.seconds(i_burst) = fread(fid,1,'uint32');
            data.microseconds(i_burst) = fread(fid, 1, 'uint32')* 1e-6;
            
            %             data_record_time = struct(  'day', day,...
%                                         'second', second, ...
%                                         'microsecond',microsecond); 
           
            %2 USO Correction Factor (ratio) 10-9 4 sl
            data.USO_correction(i_burst) = fread(fid,1,'uint32')*1e-15;
            %3 Mode ID 2 us (see table 2.3.3-2) 
            modeID_mode = fread(fid,1,'ubit6');
            modeID_SARIn_deg_case = fread(fid,1,'ubit1');
            modeID_reserved = fread(fid,1,'ubit1');
            modeID_CAL4 = fread(fid,1,'ubit1');
            modeID_Platform_Cont = fread(fid,1,'ubit2');
            modeID_reserved2 = fread(fid,1,'ubit5'); 

            modeID_CR = struct('modeID_mode', modeID_mode,...
                                    'modeID_SARIn_deg_case', modeID_SARIn_deg_case,...
                                    'modeID_reserved', modeID_reserved,...
                                    'modeID_CAL4', modeID_CAL4,...
                                    'modeID_Platform_Cont', modeID_Platform_Cont,...
                                    'modeID_reserved2',modeID_reserved2);
           switch modeID_mode
              case 2 %SAR
                   data.ProcessID(i_burst) = 58;
               case 3 %SARIN
                   data.ProcessID(i_burst) = 58;
                   if modeID_CAL4
                       
                       data.ProcessID(i_burst)= 57; %CAL4 burst
                       
                    end
           end
                   
            
            %4 Source Sequence Counter 2 us (see note 6)
            data.source_seq_count_sar_ku_fbr(i_burst) = fread(fid,1,'uint16');
            %5 Instrument Configuration 4 ul (see table 2.3.3-3)

            ins_rx_in_use = fread(fid,1,'ubit2');
            data.ins_id(i_burst) = fread(fid,1,'ubit1');
            ins_reserved = fread(fid,1,'ubit1');
            ins_bandwidth = fread(fid,1,'ubit2');
            ins_reserved2 = fread(fid,1,'ubit2');
            data.ins_tracking_mode(i_burst) = fread(fid,1,'ubit2');
            ins_ext_calibration = fread(fid,1,'ubit1');
            ins_reserved3 = fread(fid,1,'ubit1');
            data.ins_loop_stat(i_burst) = fread(fid,1,'ubit1');
            ins_loss_echo = fread(fid,1,'ubit1');
            ins_real_time_error = fread(fid,1,'ubit1');
            ins_echo_sat_error = fread(fid,1,'ubit1');
            ins_rx_band_attenuation = fread(fid,1,'ubit1');
            ins_cycle_report = fread(fid,1,'ubit1');
            ins_star_tracker_1 = fread(fid,1,'ubit1');
            ins_star_tracker_2 = fread(fid,1,'ubit1');
            ins_star_tracker_3 = fread(fid,1,'ubit1');
            ins_reserved4 = fread(fid,1,'ubit11');
            
            data.inst_id_sar_isp(i_burst)=0;
%             instrument_configuration =struct(    'ins_rx_in_use', ins_rx_in_use,...
%                                 'ins_siral_id', ins_siral_id,...
%                                 'ins_reserved', ins_reserved,...
%                                 'ins_bandwidth', ins_bandwidth,...
%                                 'ins_reserved2', ins_reserved2,...
%                                 'ins_tracking_mode',ins_tracking_mode,...                                
%                                 'ins_ext_calibration', ins_ext_calibration,...
%                                 'ins_reserved3', ins_reserved3,...
%                                 'ins_loop_stat', ins_loop_stat,...
%                                 'ins_loss_echo', ins_loss_echo,...
%                                 'ins_real_time_error',ins_real_time_error,...
%                                 'ins_echo_sat_error', ins_echo_sat_error,...
%                                 'ins_rx_band_attenuation', ins_rx_band_attenuation,...
%                                 'ins_cycle_report', ins_cycle_report,...
%                                 'ins_star_tracker_1', ins_star_tracker_1,...
%                                 'ins_star_tracker_2', ins_star_tracker_2,...
%                                 'ins_star_tracker_3', ins_star_tracker_3,...
%                                 'ins_reserved4',ins_reserved4);

            %6 Burst counter (always starts from 1 and incremented at group rate) 4 ul
            
            data.pri_sar_isp(i_burst) = pri_sar_chd;
            data.ambiguity_order_sar_isp(i_burst) = 0;
            data.burst_sar_ku(i_burst) = fread(fid,1,'uint32');
            data.burst_sar_ku_fbr(i_burst)= 1+mod(data.burst_sar_ku(i_burst),N_bursts_cycle_chd);
%             if(data.burst_sar_ku_fbr(i_burst)==0)
%                 data.burst_sar_ku_fbr(i_burst)=1;
%             end
            %7 Latitude of measurement 10-1 mcrodeg 4 sl (see note 1)
            data.lat_sar_sat(i_burst) = fread (fid,1,'int32') * 1e-7;
            
            %8 Longitude of measurement 10-1 microdeg 4 sl (see note 1)
            data.lon_sar_sat(i_burst) = fread (fid,1,'int32') * 1e-7;
            
            %9 Altitude of COG above reference ellipsoid (interpolated value)mm 4 sl
            data.alt_sar_sat(i_burst) = fread (fid,1,'int32') * 1e-3;
            
            %10 Instantaneous altitude rate derived from orbit mm/s 4 sl
             data.alt_rate_sar_sat(i_burst) = fread (fid,1,'int32') * 1e-3;
             
            
             %11 Satellite velocity vector[3](in ITRF) mm/s 3*4 sl
             data.x_vel_sat_sar(i_burst) = fread(fid,1,'int32') * 1e-3;
             
             data.y_vel_sat_sar(i_burst) = fread(fid,1,'int32') * 1e-3;
             
             data.z_vel_sat_sar(i_burst) = fread(fid,1,'int32') * 1e-3;
             
             
            %12 Real beam direction vector[3](in CRF) micros 3*4 sl
            real_beam_direction(1) = fread(fid,1,'int32')* 1e-6;
            real_beam_direction(2) = fread(fid,1,'int32')* 1e-6;
            real_beam_direction(3) = fread(fid,1,'int32')* 1e-6;
            %13 Interferometer baseline vector[3](in CRF)micros 3*4 sl
            inferometer_baseline(1) = fread(fid,1,'int32')* 1e-6;
            inferometer_baseline(2) = fread(fid,1,'int32')* 1e-6;
            inferometer_baseline(3) = fread(fid,1,'int32')* 1e-6;
            
            data.roll_sar(i_burst)  = atan((inferometer_baseline(1)/inferometer_baseline(3))); % [rad]
            data.pitch_sar(i_burst) = atan((-real_beam_direction(2)/real_beam_direction(1))); % [rad]
            data.yaw_sar(i_burst)   = atan((inferometer_baseline(2)/inferometer_baseline(3))); % [rad]
            
        
            
            %14 FBR Measurement ConfidenceData (flag word)4 u
            data.confi_block_degraded(i_burst) = fread(fid,1,'ubit1');
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
    
%                 mea_conf_data_sar_ku_fbr(i_burst)=0; %TBD
                data.mea_conf_data_sar_ku_fbr(i_burst) =struct(    'confi_block_degraded', data.confi_block_degraded(i_burst),...
                                    'confi_blank_block', confi_blank_block,...
                                    'confi_datation_degraded', confi_datation_degraded,...
                                    'confi_orbit_prop_error', confi_orbit_prop_error,...
                                    'confi_orbit_file_change', confi_orbit_file_change,...
                                    'confi_orbit_discon',confi_orbit_discon,...                                
                                    'confi_echo_sat', confi_echo_sat,...
                                    'confi_other_echo_error', confi_other_echo_error,...
                                    'confi_rx1_error_for_SARIN', confi_rx1_error_for_SARIN,...
                                    'confi_rx2_error_for_SARIN', confi_rx2_error_for_SARIN,...
                                    'confi_window_delay_inconsistency',confi_window_delay_inconsistency,...
                                    'confi_AGC_inconsistency', confi_AGC_inconsistency,...
                                    'confi_cal1_correction_miss', confi_cal1_correction_miss,...
                                    'confi_cal1_correction_from_IPF_DB', confi_cal1_correction_from_IPF_DB,...
                                    'confi_DORIS_USO_correction', confi_DORIS_USO_correction,...
                                    'confi_complex_cal1_correction_from_IPF_DB', ins_star_tracker_2,...
                                    'confi_TRK_echo_error', ins_star_tracker_3,...
                                    'confi_echo_rx1_error',ins_reserved4,...                            
                                    'confi_echo_rx2_error', confi_echo_rx2_error,...
                                    'confi_NMP_inconsistency', confi_NMP_inconsistency,...
                                    'confi_azimuth_cal_missing', confi_azimuth_cal_missing,...
                                    'confi_azimuth_cal_from_IPF_DB', confi_azimuth_cal_from_IPF_DB,...
                                    'confi_range_window_cal_function_missing',confi_range_window_cal_function_missing,...
                                    'confi_range_window_cal_function_from_IPF_DB', confi_range_window_cal_function_from_IPF_DB,...
                                    'confi_reserved', confi_reserved,...
                                    'confi_cal2_correction_missing', confi_cal2_correction_missing,...
                                    'confi_cal2_correction_from_IPF_DB', confi_cal2_correction_from_IPF_DB,...
                                    'confi_power_scaling_error_LRM_only', confi_power_scaling_error_LRM_only,...
                                    'confi_attitude_correction_missing', confi_attitude_correction_missing,...
                                    'confi_attitude_interpolation_error',confi_attitude_interpolation_error,...
                                    'confi_reserved2',confi_reserved2,...
                                    'confi_phase_perturbation',confi_phase_perturbation);

%                 time_orbit_group(j)=struct(    'data_record_time', data_record_time,...
%                                         'USO_correction', USO_correction,...
%                                         'modeID', modeID,...
%                                         'ssc', ssc,...
%                                         'instrument_configuration', instrument_configuration,...
%                                         'burst_counter',burst_counter,...
%                                         'latitude',latitude,...
%                                         'longitude',longitude,...
%                                         'altitude_COG',altitude_COG,...
%                                         'instantaneous_altitude_rate',instantaneous_altitude_rate,...
%                                         'satellite_velocity',satellite_velocity,...
%                                         'real_beam_direction',real_beam_direction,...
%                                         'inferometer_baseline',inferometer_baseline,...
%                                         'FBR_measurement_confidence_flag',FBR_measurement_confidence_flag);

    
        end
        
    end
    
    %-----------------------------%
    %-- read Measurements Group --%
    %-----------------------------%
          
    
    for j=1:20
        
        i_burst= (record_num-1)*20+j;
                    
		%15 Window Delay (2way) uncorrected for instrument delays 10-12 s 8 sll
        data.win_delay_sar_ku(i_burst) = fread(fid,1,'int64') * 1e-12;
		%16 Ho Initial Height Word 48.8 ps 4 sl (see note 2)
        data.h0_comp_sar_isp(i_burst)  = fread(fid,1,'int32') * 1e-12 * 48.8;
        %17 HPR Height Rate 3.05 ps/rc 4 sl (see note 2)
        data.cor2_comp_sar_isp(i_burst) = fread(fid,1,'int32') * 1e-12* 3.05;
		%18 LAI 12.5 ns 4 sl (see note 2)
        data.h0_sar_isp(i_burst) = fread(fid,1,'int32') * 1e-9 * 12.5;
		%19 FAI 12.5/256 ns 4 sl (see note 2)
        data.cor2_sar_isp(i_burst) = fread(fid,1,'int32') * 1e-9 * 12.5 / 256;
		%20 AGC_1 (not corrected) dB/100 4 sl (see note 3)
        data.ATT1_science(i_burst) = -1.0 * fread(fid,1,'int32')/100;
		%21 AGC_2 (not corrected) dB/100 4 sl (see note 3)
        data.ATT2_science(i_burst) = -1.0 * fread(fid,1,'int32')/100;
        data.att_sar_ku_isp(i_burst) = 62-(data.ATT1_science(i_burst));
            
        
        
		%22 Total Fixed Gain Rx 1 dB/100 4 sl (see note 3)
        data.tot_fixed_gain_1(i_burst) = fread(fid,1,'int32')/100;
		%23 Total Fixed Gain Rx 2 dB/100 4 sl (see note 3)
        data.tot_fixed_gain_2(i_burst) = fread(fid,1,'int32')/100;
		%24 Transmit Power Micro-Watts 4 sl
        data.transmit_power(i_burst) = fread(fid,1,'int32')* 1e-6;
		%25 Doppler range correction (Radial component)mm 4 sl
        data.doppler_range_correction(i_burst) = fread(fid,1,'int32') * 1e-3;
		%26 Instrument Range Correction tx-rx antenna mm 4 sl
        data.instrument_range_correction_tx_rx(i_burst) = fread(fid,1,'int32') * 1e-3;
		%27 Instrument Range Correction mm 4 sl
        data.instrument_range_correction(i_burst) = fread(fid,1,'int32') * 1e-3;
		%28 Instrument Sigma 0 correction, tx-rx antenna dB/100 4 sl (see note
        data.instrument_sigma0_correction_tx_rx(i_burst) = fread(fid,1,'int32')/100;
		%29 Instrument Sigma 0 correction rx only antenna dB/100 4 sl (see note
        data.instrument_sigma0_correction_rx(i_burst) = fread(fid,1,'int32')/100;
		%30 Internal Phase Correction Microradians 4 sl
        data.internal_phase_correction(i_burst) = fread(fid,1,'int32') * 1e-6;
		%31 External Phase Correction Microradians 4 sl
        data.external_phase_correction(i_burst) = fread(fid,1,'int32')  * 1e-6;
		%32 Noise power measurement dB/100 4 sl (see note
        data.noise_power(i_burst) = fread(fid,1,'int32') /100;
		%33 Phase Slope Correction Microradians 4 sl (see note
        data.phase_slope_correction(i_burst) = fread(fid,1,'int32') * 1e-6;
        %34 spares 4*1 uc
        fread(fid,4,'uchar');
        
        data.gain_corr_instr_sar(i_burst)= data.tot_fixed_gain_1(i_burst)-(data.ATT1_science(i_burst)+data.ATT2_science(i_burst))+data.instrument_sigma0_correction_tx_rx(i_burst);

%         measurements_group(j)=struct(   'window_delay',window_delay,...
%                                         'ho_initial_height',ho_initial_height,...
%                                         'HPR_hight_rate',HPR_hight_rate,...
%                                         'LAI',LAI,...
%                                         'FAI',FAI,...
%                                         'AGC_1',AGC_1,...
%                                         'AGC_2',AGC_2,...
%                                         'tot_fixed_gain_1',tot_fixed_gain_1,...
%                                         'tot_fixed_gain_2',tot_fixed_gain_2,...
%                                         'transmit_power',transmit_power,...
%                                         'doppler_range_correction',doppler_range_correction,...
%                                         'instrument_range_correction_tx_rx',instrument_range_correction_tx_rx,...
%                                         'instrument_range_correction',instrument_range_correction,...
%                                         'instrument_sigma0_correction_tx_rx',instrument_sigma0_correction_tx_rx,...
%                                         'instrument_sigma0_correction_rx',instrument_sigma0_correction_rx,...
%                                         'internal_phase_correction',internal_phase_correction,...
%                                         'external_phase_correction',external_phase_correction,...
%                                         'noise_power',noise_power);
                                            
%       disp(['Measurements Group ' num2str(j) '/' num2str(i)]);
% 		disp(['window_delay ' num2str(window_delay)]);
% 		disp(['ho_initial_height ' num2str(ho_initial_height)]);
% 		disp(['HPR_hight_rate ' num2str(HPR_hight_rate)]);
% 		disp(['LAI ' num2str(LAI)]);
% 		disp(['FAI ' num2str(FAI)]);
% 		disp(['AGC_1 ' num2str(AGC_1)]);
% 		disp(['AGC_2 ' num2str(AGC_2)]);
% 		disp(['tot_fixed_gain_1 ' num2str(tot_fixed_gain_1)]);
% 		disp(['tot_fixed_gain_2 ' num2str(tot_fixed_gain_2)]);
% 		disp(['transmit_power ' num2str(transmit_power)]);
% 		disp(['doppler_range_correction ' num2str(doppler_range_correction)]);
% 		disp(['instrument_range_correction_tx_rx ' num2str(instrument_range_correction_tx_rx)]);
% 		disp(['instrument_range_correction ' num2str(instrument_range_correction)]);
% 		disp(['instrument_sigma0_correction_tx_rx ' num2str(instrument_sigma0_correction_tx_rx)]);
% 		disp(['instrument_sigma0_correction_rx ' num2str(instrument_sigma0_correction_rx)]);
% 		disp(['internal_phase_correction ' num2str(internal_phase_correction)]);
% 		disp(['external_phase_correction ' num2str(external_phase_correction)]);
% 		disp(['noise_power ' num2str(noise_power)]);
        
    end
    
    if(end_of_file==0)
        %-----------------------------%
        %--  read Corrections Group --%
        %-----------------------------%
        %35 Dry Tropospheric Correction mm 4 sl
        data.dry_tropo_correction(record_num) = fread(fid,1,'int32')  * 1e-3;
        data.dry_tropo_correction_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.dry_tropo_correction(record_num);
        %36 Wet Tropospheric Correction mm 4 sl
        data.wet_tropo_correction(record_num) = fread(fid,1,'int32') * 1e-3;
        data.wet_tropo_correction_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.wet_tropo_correction(record_num);
        %37 Inverse Barometric Correction mm 4 sl
        data.inverse_baro_correction(record_num) = fread(fid,1,'int32') * 1e-3;
        data.inverse_baro_correction_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.inverse_baro_correction(record_num);
        %38 Dynamic Atmospheric Correction mm 4 sl
        data.Dynamic_atmospheric_correction(record_num) = fread(fid,1,'int32') * 1e-3;        
        data.Dynamic_atmospheric_correction_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.Dynamic_atmospheric_correction(record_num);
        %39 GIM Ionospheric Correction mm 4 sl
        data.GIM_iono_correction(record_num) = fread(fid,1,'int32') * 1e-3;
        data.GIM_iono_correction_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.GIM_iono_correction(record_num);
        %40 Model Ionospheric Correction mm 4 sl
        data.model_iono_correction(record_num) = fread(fid,1,'int32') * 1e-3;
        data.model_iono_correction_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.model_iono_correction(record_num);
        %41 Ocean Equilibrium Tide mm 4 sl
        data.ocean_equilibrium_tide(record_num) = fread(fid,1,'int32') * 1e-3;
        data.ocean_equilibrium_tide_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.ocean_equilibrium_tide(record_num);
        %42 Long Period Tide Height 4 sl
        data.long_period_tide_height(record_num) = fread(fid,1,'int32') * 1e-3;
        data.long_period_tide_height_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.long_period_tide_height(record_num);
        %43 Ocean Loading Tide mm 4 sl
        data.ocean_loading_tide(record_num) = fread(fid,1,'int32') * 1e-3;
        data.ocean_loading_tide_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.ocean_loading_tide(record_num);
        %44 Solid Earth Tide mm 4 sl
        data.solid_earth_tide(record_num) = fread(fid,1,'int32') * 1e-3;
        data.solid_earth_tide_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.solid_earth_tide(record_num);
        %45 Geocentric Polar Tide mm 4 sl
        data.geocentric_polar_tide(record_num) = fread(fid,1,'int32') * 1e-3;
        data.geocentric_polar_tide_bursts((record_num-1)*20+1:(record_num-1)*20+20)=data.geocentric_polar_tide(record_num);
        %46 Surface type flag 4 ul
        data.surface_type_flag(record_num) = fread(fid,1,'uint32');
        data.surface_type_flag_bursts((record_num-1)*20+1:(record_num-1)*20+20)= data.surface_type_flag(record_num);
        %47 spares 4*1 uc
        fread(fid,4,'uint8');
        %48 Correction status flags 4 ul (see table 2.3.3-5)
        %correction_status_flags = fread(fid,1,'uint32');

        correction_status_flags(1) = fread(fid,1,'ubit1');
        correction_status_flags(2) = fread(fid,1,'ubit1');
        correction_status_flags(3) = fread(fid,1,'ubit1');
        correction_status_flags(4) = fread(fid,1,'ubit1');
        correction_status_flags(5) = fread(fid,1,'ubit1');
        correction_status_flags(6) = fread(fid,1,'ubit1');
        correction_status_flags(7) = fread(fid,1,'ubit1');
        correction_status_flags(8) = fread(fid,1,'ubit1');
        correction_status_flags(9) = fread(fid,1,'ubit1');
        correction_status_flags(10) = fread(fid,1,'ubit1');
        correction_status_flags(11) = fread(fid,1,'ubit1');
        correction_status_flags(12) = fread(fid,1,'ubit1');
        correction_status_flags(13) = fread(fid,1,'ubit20');
                
        
        %49 Correction error flags 4 ul (see table 2.3.3-6)
        %correction_error_flags = fread(fid,1,'uint32');
        
        correction_error_flags(1) = fread(fid,1,'ubit1');
        correction_error_flags(2) = fread(fid,1,'ubit1');
        correction_error_flags(3) = fread(fid,1,'ubit1');
        correction_error_flags(4) = fread(fid,1,'ubit1');
        correction_error_flags(5) = fread(fid,1,'ubit1');
        correction_error_flags(6) = fread(fid,1,'ubit1');
        correction_error_flags(7) = fread(fid,1,'ubit1');
        correction_error_flags(8) = fread(fid,1,'ubit1');
        correction_error_flags(9) = fread(fid,1,'ubit1');
        correction_error_flags(10) = fread(fid,1,'ubit1');
        correction_error_flags(11) = fread(fid,1,'ubit1');
        correction_error_flags(12) = fread(fid,1,'ubit1');
        correction_error_flags(13) = fread(fid,1,'ubit20');
        %50 Spare 4*1 uc
        fread(fid,4,'uint8');


%         corrections_group = struct( 'dry_tropo_correction',dry_tropo_correction,...
%                                     'wet_tropo_correction',wet_tropo_correction,...
%                                     'inverse_baro_correction',inverse_baro_correction,...
%                                     'Dynamic_atmospheric_correction',Dynamic_atmospheric_correction,...
%                                     'GIM_iono_correction',GIM_iono_correction,...
%                                     'model_iono_correction',model_iono_correction,...
%                                     'ocean_equilibrium_tide',ocean_equilibrium_tide,...
%                                     'long_period_tide_height',long_period_tide_height,...
%                                     'ocean_loading_tide',ocean_loading_tide,...
%                                     'solid_earth_tide',solid_earth_tide,...
%                                     'geocentric_polar_tide',geocentric_polar_tide,...
%                                     'surface_type_flag',surface_type_flag,...
%                                     'correction_status_flags',correction_status_flags,...
%                                     'correction_error_flags',correction_error_flags);

        %--------------------%
        %-- waveform group --%
        %--------------------%

  
%   fread(fid,2621520,'int8'); 
%   SAR FBR waveform group 32760  
%   SARIn FBR waveforms group 2621520 
    for j=1:20
        
        i_burst= (record_num-1)*20+j;
        progressbar(i_burst/(num_dsr*20+1),[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
        % 54--55 Complex Echo Waveform [64,512,2] 512*2 sc 
%             for iSample=1:128
%                 averaged_power_echo_waveform(iSample) = fread(fid,1,'int16');
%             end
            wfm_iq_sar_ku_fbr_aux_1(:) = fread (fid, N_ku_pulses_burst_chd * N_samples_sar_chd*2 ,'int8');
%             wfm_iq_sar_ku_fbr_aux_2(:) = fread (fid, N_ku_pulses_burst_chd * N_samples_sar_chd*2 ,'int8');
            
            wfm_iq_sar_ku_fbr_i_1 = wfm_iq_sar_ku_fbr_aux_1(1:2:N_ku_pulses_burst_chd * N_samples_sar_chd*2);
            wfm_iq_sar_ku_fbr_q_1 = wfm_iq_sar_ku_fbr_aux_1(2:2:N_ku_pulses_burst_chd * N_samples_sar_chd*2);
%             wfm_iq_sar_ku_fbr_i_2 = wfm_iq_sar_ku_fbr_aux_2(1:2:N_ku_pulses_burst_chd * N_samples_sar_chd*2);
%             wfm_iq_sar_ku_fbr_q_2 = wfm_iq_sar_ku_fbr_aux_2(2:2:N_ku_pulses_burst_chd * N_samples_sar_chd*2);
            
            for i_pulse = 1:N_ku_pulses_burst_chd
                wfm_iq_sar_ku_fbr_i_11(i_burst,i_pulse,:) = wfm_iq_sar_ku_fbr_i_1((i_pulse-1)*N_samples_sar_chd+1:i_pulse*N_samples_sar_chd)./(10.^(data.gain_corr_instr_sar(i_burst)./20));
                wfm_iq_sar_ku_fbr_q_11(i_burst,i_pulse,:) = wfm_iq_sar_ku_fbr_q_1((i_pulse-1)*N_samples_sar_chd+1:i_pulse*N_samples_sar_chd)./(10.^(data.gain_corr_instr_sar(i_burst)./20));
                %                 wfm_iq_sar_ku_fbr_i_22 = wfm_iq_sar_ku_fbr_i_2((i_pulse-1)*N_samples_sar_chd+1:i_pulse*N_samples_sar_chd);
                %                 wfm_iq_sar_ku_fbr_q_22 = wfm_iq_sar_ku_fbr_q_2((i_pulse-1)*N_samples_sar_chd+1:i_pulse*N_samples_sar_chd);
            end           
            
            % 52 Echo Scale Factor (to scale echo to watts) - 4 sl
            data.nimp_sar_isp(i_burst) = fread(fid,1,'uint16');
%
            fread(fid,1,'uint16');
    end

    end
end
clear wfm_iq_sar_ku_fbr_i_1 wfm_iq_sar_ku_fbr_q_1
% data.wfm_sar_reversed = wfm_iq_sar_ku_fbr_i_2 + 1i * wfm_iq_sar_ku_fbr_q_2;

data.USO_correction = data.win_delay_sar_ku.*data.USO_correction; %output file in seconds
data.win_delay_sar_ku = data.win_delay_sar_ku + data.USO_correction + (data.instrument_range_correction_tx_rx )/(c_cst/2);

% figure; imagesc(abs(fftshift(fft(squeeze(wfm_iq_sar_ku_fbr_i_2(1,:,:)-1i*wfm_iq_sar_ku_fbr_q_2(1,:,:))'),1)))
% figure; imagesc(abs(fftshift(fft(squeeze(wfm_iq_sar_ku_fbr(1,:,:))'),2)))

data.wfm_cal_gain_corrected = wfm_iq_sar_ku_fbr_i_11 + 1i * wfm_iq_sar_ku_fbr_q_11;
%data.wfm_cal_gain_corrected_2 = wfm_iq_sar_ku_fbr_i_22 + 1i * wfm_iq_sar_ku_fbr_q_22;
clear wfm_iq_sar_ku_fbr_i_2 wfm_iq_sar_ku_fbr_q_2 
fclose(fid);

data.N_total_bursts_sar_ku = length(data.lat_sar_sat);
data.time_sar_ku = data.days + data.seconds + data.microseconds;

%Modified by EM 30.03.2016: Include CAL1-pulse-to-pulse and CAL2
%calibrations
for i_burst=data.N_total_bursts_sar_ku
    %------------------- CAL1-pulse-to-pulse-----------------------
    %operating using matrices instead of loops
    if CAL1pp_active
        %--------- Apply CAL1-pulse-to-pulse calibration ----------
        aux=squeeze(data.wfm_cal_gain_corrected(i_burst,:,:));
        aux=aux.*((10.^ (burst_power_array_cor_cal1_sar_rep.'/20))*ones(1,N_samples_sar_chd)).*...
            exp (1i * (burst_phase_array_cor_cal1_sar_rep.'*ones(1,N_samples_sar_chd)));
        data.wfm_cal_gain_corrected(i_burst,:,:) = aux;
        clear aux;
%         aux=squeeze(data.wfm_cal_gain_corrected_2(i_burst,:,:));
%         aux=aux.*((10 ^ (burst_power_array_cor_cal1_sar_rep.'/20))*ones(1,N_samples_sar_chd)).*...
%             exp (1i * (burst_phase_array_cor_cal1_sar_rep.'*ones(1,N_samples_sar_chd)));
%         data.wfm_cal_gain_corrected_2(i_burst,:,:) = aux;
%         clear aux;
                    
    end
    %------------------ CAL2---------------------------------------
    %operating using matrices instead of loops
    if CAL2_active
        aux=fft(squeeze(data.wfm_cal_gain_corrected(i_burst,:,:)),N_samples_sar_chd,2); % fft along range samples
        data.wfm_cal_gain_corrected(i_burst,:,:)=ifft(aux.*(ones(N_ku_pulses_burst_chd,1)*wfm_cal2_science_sar_rep),N_samples_sar_chd,2);
        clear aux;
%         aux=fft(squeeze(data.wfm_cal_gain_corrected(i_burst,:,:)),N_samples_sar_chd,2); % fft along range samples
%         data.wfm_cal_gain_corrected_2(i_burst,:,:)=ifft(aux.*(ones(N_ku_pulses_burst_chd,1)*wfm_cal2_science_sar_rep),N_samples_sar_chd,2);
%         clear aux;
    end
end

p = lla2ecef([data.lat_sar_sat.',data.lon_sar_sat.',data.alt_sar_sat.'],flat_coeff_cst,semi_major_axis_cst);
data.x_sar_sat= p(:,1).';
data.y_sar_sat = p(:,2).';
data.z_sar_sat = p(:,3).';

alt_surf_out = data.alt_sar_sat - data.win_delay_sar_ku * c_cst/2;

% geod2cart(SURF)
p = lla2ecef([data.lat_sar_sat',data.lon_sar_sat',alt_surf_out'],flat_coeff_cst,semi_major_axis_cst);
x_surf_geoloc = p(:,1);
y_surf_geoloc = p(:,2);
z_surf_geoloc = p(:,3);

% data=filter_data_err_FBR_CrypSat_SAR(data);    
    

[~,data.doppler_ang_sar_sat] = compute_height_rate(data.N_total_bursts_sar_ku, data.x_vel_sat_sar, data.y_vel_sat_sar, data.z_vel_sat_sar,data.x_sar_sat ,data.y_sar_sat ,data.z_sar_sat ,x_surf_geoloc,y_surf_geoloc,z_surf_geoloc);
data.T0_sar     = T0_chd.* ones(1,data.N_total_bursts_sar_ku);
pri_sar_nom = pri_sar_chd * ones(1,data.N_total_bursts_sar_ku);
data.pri_sar = pri_sar_nom;
prf_sar_nom = prf_sar_chd * ones(1,data.N_total_bursts_sar_ku); 
bri_nom = bri_chd * ones(1,data.N_total_bursts_sar_ku); 
brf_nom = brf_chd * ones(1,data.N_total_bursts_sar_ku); 



temps = toc(t1);
minuts = floor(temps/60);
segs = temps - minuts*60;
disp(['Han passat ',num2str(minuts),' minuts i ',num2str(segs),' segons']);

end