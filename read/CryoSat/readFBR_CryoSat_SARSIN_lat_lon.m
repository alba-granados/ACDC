% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% 
%
% ---------------------------------------------------------
% READFBR: function that reads the FBR data set from input filename,
% assuming SAR or SARIN data records
%
% Calling
%   fbr_ds = readFBR( filename, headers )
%
% Inputs
%   filename: input SAR SARIN FBR file
%   headers:   if false --> without header
%              if true --> with header
%
% Output
%   fbr_ds           : data contained in the file 
%
% ----------------------------------------------------------
% 
% Author:   Albert Garcia / isardSAT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [lat_sar_sat,lon_sar_sat, alt_sar_sat, win_delay_sar_ku, time_sar_ku] = readFBR_CryoSat_SARSIN_lat_lon(filename_FBR)


global sec_in_day_cst N_bursts_cycle_chd
global pri_sar_chd 

sph = readSPH(filename_FBR);
% FBR data set must be the first one
fbr_dsd = sph.dsds(1);
ds_offset = fbr_dsd.ds_offset;
ds_size = fbr_dsd.ds_size;
 num_dsr = fbr_dsd.num_dsr;
dsr_size = fbr_dsd.dsr_size;

lat_sar_sat = zeros(1,num_dsr*20); 
lon_sar_sat = zeros(1,num_dsr*20);
alt_sar_sat = zeros(1,num_dsr*20);
win_delay_sar_ku = zeros(1,num_dsr*20);

%% From Product Specs
% Fields mapping size           1   2 	3   4   5   6   7   8   9   10  11      12      13      14
% time_orbit_group_record_size =  12+ 4+  2+  2+  4+  4+  4+  4+  4+  4+  3*4+    3*4+    3*4+    4;
% Fields mapping size           15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31	32 	33	34
measurements_group_record_size = 8+ 4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4;
% Fields mapping size           35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  
corrections_group_record_size  = 4+ 4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4;

if(strfind(fbr_dsd.ds_name,'SIR_FBR_SARIN'))
    % Fields mapping size           54          55          56  57
    waveform_group_record_size     = 64*512*2+  64*512*2+   2+  2;%2621520/20
elseif(strfind(fbr_dsd.ds_name,'SIR_FBR_SAR'))
    % Fields mapping size           51          52  53 
    waveform_group_record_size     = 64*128*2+  2+  2; %327760/20

end


% read the dsrs
fid = fopen(filename_FBR,'r','b');
% fseek(fid,-num_dsr*dsr_size,'eof'); % last records of the product
fseek(fid,ds_offset,'bof'); % start reading from the first record




% modes (by modeID_record, 32bits, TABLE02)

% progressbar('Records','Data Blocks')

%% READ FBR FILE 
for record_num=1:num_dsr
    %--------------------------------%
    %--  read Time and Orbit Group --%
    %--------------------------------%
    
    for j=1:20

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
            USO_correction = fread(fid,1,'uint32')*1e-15;
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
            ins_siral_id = fread(fid,1,'ubit1');
            ins_reserved = fread(fid,1,'ubit1');
            ins_bandwidth = fread(fid,1,'ubit2');
            ins_reserved2 = fread(fid,1,'ubit2');
            ins_tracking_mode = fread(fid,1,'ubit2');
            ins_ext_calibration = fread(fid,1,'ubit1');
            ins_reserved3 = fread(fid,1,'ubit1');
            ins_loop_stat = fread(fid,1,'ubit1');
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

            %6 Burst counter (always starts from 1 and incremented at group rate) 4 ul
            
            data.pri_sar_isp(i_burst) = pri_sar_chd;
            data.ambiguity_order_sar_isp(i_burst) = 0;
            data.burst_sar_ku(i_burst) = fread(fid,1,'uint32');
            data.burst_sar_ku_fbr(i_burst)= mod(data.burst_sar_ku(i_burst),N_bursts_cycle_chd);
            if(data.burst_sar_ku_fbr(i_burst)==0)
                data.burst_sar_ku_fbr(i_burst)=1;
            end
            %7 Latitude of measurement 10-1 mcrodeg 4 sl (see note 1)
            lat_sar_sat(i_burst) = fread (fid,1,'int32') * 1e-7;
            
            %8 Longitude of measurement 10-1 microdeg 4 sl (see note 1)
            lon_sar_sat(i_burst) = fread (fid,1,'int32') * 1e-7;
            
            %9 Altitude of COG above reference ellipsoid (interpolated value)mm 4 sl
            alt_sar_sat(i_burst) = fread (fid,1,'int32') * 1e-3;
            
            %10 Instantaneous altitude rate derived from orbit mm/s 4 sl
             data.alt_rate_sar_sat(i_burst) = fread (fid,1,'int32') * 1e-3;
             
            
             %11 Satellite velocity vector[3](in ITRF) mm/s 3*4 sl
             data.x_vel_sat_sar(i_burst) = fread(fid,1,'int32') * 1e-3;
             
             data.y_vel_sat_sar(i_burst) = fread(fid,1,'int32') * 1e-3;
             
             data.z_vel_sat_sar(i_burst) = fread(fid,1,'int32') * 1e-3;
             
             
            %12 Real beam direction vector[3](in CRF) micros 3*4 sl
            real_beam_direction(1) = fread(fid,1,'int32');
            real_beam_direction(2) = fread(fid,1,'int32');
            real_beam_direction(3) = fread(fid,1,'int32');
            %13 Interferometer baseline vector[3](in CRF)micros 3*4 sl
            inferometer_baseline(1) = fread(fid,1,'int32');
            inferometer_baseline(2) = fread(fid,1,'int32');
            inferometer_baseline(3) = fread(fid,1,'int32');
            
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
%                 data.mea_conf_data_sar_ku_fbr(i_burst) =struct(    'confi_block_degraded', data.confi_block_degraded(i_burst),...
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
    
    %-----------------------------%
    %-- read Measurements Group --%
    %-----------------------------%
	for j=1:20
        i_burst= (record_num-1)*20+j;
		%15 Window Delay (2way) uncorrected for instrument delays 10-12 s 8 sll
        win_delay_sar_ku(i_burst) = fread(fid,1,'int64') * 1e-12;
        fseek(fid,(measurements_group_record_size-8),'cof');
    end
        
    %-----------------------------%
    %--  read Corrections Group --%
    %-----------------------------%
    fseek(fid,corrections_group_record_size,'cof'); 

    %--------------------%
    %-- waveform group --%
    %--------------------%

    fseek(fid,waveform_group_record_size*20,'cof'); 


end

fclose(fid);

time_sar_ku = data.days + data.seconds + data.microseconds;


data.N_total_bursts_sar_ku = length(data.days);



end