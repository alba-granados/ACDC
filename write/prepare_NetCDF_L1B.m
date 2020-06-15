%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the CODING & PACKING 
% algorithm as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v6c_20140604
% and JC-ID-ESA-GP-0059_issue20_140704 (IODD)
%
% ---------------------------------------------------------
% Objective: Pack variables and write the NETCDF
% 
% INPUTs : Workspace
% OUTPUTs: TM Structure as defined on isardSAT_JasonCS_DPM
%
% ----------------------------------------------------------
% Author:    Gorka Moyano  / isardSAT
%            Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/07/2014)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils

t6 = tic;
  no_amb = 0;
% no_amb = 1; Waveforms without ambiguities in the product
% no_amb = 0; Waveforms with ambiguities in the product


if(no_amb)
%     ambig_text='NoAmb'; % Ambiguities filtered
    ambig_text='2'; % Ambiguities filtered
else
%     ambig_text='Amb'; % Ambiguities not removed
    ambig_text='1'; % Ambiguities not removed
end

if strcmp(TEST,'S6_C1A_SAR')
    filename = ['S6_P4_SIM_RAW_L1B__20190119T064000_20190119T064019_T01' ambig_text '.nc'];
elseif strcmp(TEST,'S6_C1A_RMC')
    filename = ['S6_P4_SIM_RMC_L1B__20190119T064000_20190119T064019_T01' ambig_text '.nc'];    
elseif strcmp(TEST,'S6_C1B_SAR')
    filename = ['S6_P4_SIM_RAW_L1B__20210929T064000_20210929T064019_T02' ambig_text '.nc'];   
elseif strcmp(TEST,'S6_C1B_RMC')
    filename = ['S6_P4_SIM_RMC_L1B__20210929T064000_20210929T064019_T02' ambig_text '.nc'];    
elseif strcmp(TEST,'S6_C1C_SAR')
    filename = ['S6_P4_SIM_RAW_L1B__20210929T064000_20210929T064019_T03' ambig_text '.nc'];      
elseif strcmp(TEST,'S6_C1C_RMC')
    filename = ['S6_P4_SIM_RMC_L1B__20210929T064000_20210929T064019_T03' ambig_text '.nc'];     
end
% filename = [TEST '_L1B_' ambig_text '.nc'];

dimensions_key = 'Dimensions';
format_key = 'Format';
data_type_key = 'DataType';

netcdf_v4_format = 'netcdf4';

ku_rec_dimension = 'Ku_rec';
nl_dimension = 'Nl';
nl_dimension_size = 488;
ns_dimension = 'Ns';
space_3D_dimension = 'space_3D';
space_3D_dimension_size = 3;

long_name_att = 'long_name';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';

int8_type = 'int8';
uint8_type = 'uint8';
int16_type = 'int16';
uint16_type = 'uint16';
int32_type = 'int32';
uint32_type = 'uint32';
uint64_type = 'uint64';
float_type = 'single';

day_units = 'day';
seconds_units = 'seconds';
number_units = '#';
degrees_units = 'degrees';
meters_units = 'meters';
meters_per_second_units = 'm/s';
dB_units = 'dB';
Hz_units = 'Hz';
T0d64_units = 'T0/64';
T0d16d64_units = 'T0/16/64';
W_per_count_units = 'Watt/#';
rad_units = 'rad';

%% CODING L1B

%Compute nadir burst for a given surface
[~,beam_index_nadir]=min(abs(-pi/2+beam_ang_surf'));

burst_index_nadir = zeros(N_total_surf_loc);

for i_surf = 1:N_total_surf_loc
    burst_index_nadir(i_surf)= burst_index(i_surf,beam_index_nadir(i_surf));
end
int_delay_cor_cal1_surf = zeros(N_total_surf_loc,1);
z_cog_ant_corr_surf = zeros(N_total_surf_loc,1);
att_conv_sar_ku_surf = zeros(N_total_surf_loc,1);
power_var_cal1_sar_ku_surf = zeros(N_total_surf_loc,1);
gain_corr_instr_sar_surf = zeros(N_total_surf_loc,1);
T0_sar_surf_nadir = zeros(N_total_surf_loc,1);
h0_sar_isp_surf = zeros(N_total_surf_loc,1);
cor2_sar_isp_surf = zeros(N_total_surf_loc,1);
mean_cai_fai_sar_surf = zeros(N_total_surf_loc,1);
pri_sar_isp_surf = zeros(N_total_surf_loc,1);
l1_instrument_configuration_ku_surf = zeros(N_total_surf_loc,1);
source_seq_count_sar_isp_surf = zeros(N_total_surf_loc,1);
source_seq_count_sar_isp= 0:(N_total_bursts_sar_ku_isp-1);
for i_surf = 1:N_total_surf_loc
    int_delay_cor_cal1_surf(i_surf)     = int_delay_cor_cal1(burst_index_nadir(i_surf));
    z_cog_ant_corr_surf(i_surf)         = z_cog_ant_corr(burst_index_nadir(i_surf));
    att_conv_sar_ku_surf(i_surf)        = att_conv_sar_ku(burst_index_nadir(i_surf));
    power_var_cal1_sar_ku_surf(i_surf)  = power_var_cal1_sar_ku(burst_index_nadir(i_surf));
    gain_corr_instr_sar_surf(i_surf)    = gain_corr_instr_sar(burst_index_nadir(i_surf));
    T0_sar_surf_nadir(i_surf)           = T0_sar_surf(i_surf,beam_index_nadir(i_surf));
    h0_sar_isp_surf(i_surf)             = h0_sar_isp(burst_index_nadir(i_surf));
    cor2_sar_isp_surf(i_surf)           = cor2_sar_isp(burst_index_nadir(i_surf));
    mean_cai_fai_sar_surf(i_surf)       = mean_cai_fai_sar(burst_index_nadir(i_surf));
    pri_sar_isp_surf(i_surf)            = pri_sar_isp(burst_index_nadir(i_surf));
    source_seq_count_sar_isp_surf(i_surf) = source_seq_count_sar_isp(burst_index_nadir(i_surf));
    
end 
%----------A. Time variables ----------------------------------------------
if (Process_ID2 == 58) % TM_ECHO_SAR
    l1_mode_id_ku = uint8(zeros(N_total_surf_loc,1)+5); % SAR Raw L1B Ku
elseif (Process_ID2 == 59) % TM_ECHO_RMC
    l1_mode_id_ku = uint8(zeros(N_total_surf_loc,1)+8); % SAR RMC L1B Ku
else % unknown
    l1_mode_id_ku = uint8(zeros(N_total_surf_loc,1)+0);
end
time_day_ku = uint32(time_surf / sec_in_day_cst);
time_seconds_ku = uint64((time_surf - double(time_day_ku) * sec_in_day_cst) * 1e6);
tm_source_sequence_counter_ku = uint16(source_seq_count_sar_isp_surf);
l1b_record_counter_ku = uint16(0:(length(win_delay_surf)-1));

%----------B. Orbit and attitude variables ------------------------------
latitude_ku = int32(lat_sat * 1e7);
for i_lon = 1:size(lon_sat,1)
    if (lon_sat(i_lon) > 180)
        longitude_ku = uint32((lon_sat +360)* 1e7); 
    else
        longitude_ku = uint32((lon_sat )* 1e7); 
    end
end
com_altitude_ku = uint32((alt_sat - 1.3e6) * 1e4);
com_altitude_rate_ku = int32(alt_rate_sat * 1e4);
com_velocity_vector_ku = int32([x_vel_sat * 1e4; y_vel_sat * 1e4; z_vel_sat * 1e4]);
satellite_mispointing_ku = int32([pitch_surf * 180/pi_cst * 1e7; roll_surf * 180/pi_cst * 1e7; yaw_surf * 180/pi_cst * 1e7]);
pitch_bias = 0; % to CHD file
roll_bias = 0; % to CHD file
yaw_bias = 0; % to CHD file
mispointing_bias_ku = int32([pitch_bias * 180/pi_cst * 1e7, roll_bias * 180/pi_cst * 1e7, yaw_bias * 180/pi_cst * 1e7]);

%----------C. Configuration and quality variables -------------------------  
l1_instrument_configuration_ku = zeros(N_total_bursts_sar_ku_isp,1);
for i_burst=1:N_total_bursts_sar_ku_isp
    bit0    = dec2bin(inst_id_sar_isp(i_burst));
    bit15   = dec2bin(tm_mode_id_sar_isp(i_burst),4);
    trk_bin = dec2bin(trk_config_sar_isp(i_burst),8);
    bit68   = trk_bin((1:3));
    bit9    = dec2bin(loss_track_criterion_isp(i_burst));
    bit10   = dec2bin(nav_bulletin_status_isp(i_burst));
    bit11to15= '0000';
     l1_instrument_configuration_ku(i_burst,1) = bin2dec([bit10 bit9 bit68 bit15 bit0]);
%     l1_instrument_configuration_ku(i_burst,1) = bin2dec([bit11to15 bit10 bit9 bit68 bit15 bit0]);
%    l1_instrument_configuration_ku(i_burst,1) = bin2dec([bit0 bit15 bit68 bit9 bit10 bit11to15]);
end

for i_surf = 1:N_total_surf_loc
    l1_instrument_configuration_ku_surf(i_surf)     = l1_instrument_configuration_ku(burst_index_nadir(i_surf));
end

l1b_mcd_ku = zeros(N_total_surf_loc,1); % TBD

%----------D. Altimeter range variables -----------------------------------+
altimeter_range_calibrated_ku = uint32(((win_delay_surf) * c_cst/2 - 1.3e6) * 1e4);
% altimeter_range_calibrated_ku = uint32(((zeros(1,N_total_surf_loc)+win_delay_surf(1)) * c_cst/2 - 1.3e6) * 1e4);
range_corr_internal_delay_ku = int16(int_delay_cor_cal1_surf * c_cst/2 * 1e4);
range_corr_com_ku = int16(-(z_cog_ant_corr_surf.') * 1e4);

%----------E. Altimeter power variables -----------------------------------
attenuator_calibrated_ku = int16(att_conv_sar_ku_surf * 1e2);
altimeter_power_drift_ku  = int16(power_var_cal1_sar_ku_surf * 1e2);
power_corr_digital_processing_ku = int16(onboard_proc_sar * 1e2);
power_scaling_to_antenna_ku = int16(gain_corr_instr_sar_surf * 1e2);

%----------F. Altimeter engineering variables -----------------------------
altimeter_clock_ku = int32(1./T0_sar_surf_nadir - 3.95e8)* 1e9;
tm_h0_ku = uint32(h0_sar_isp_surf);
tm_cor2_ku = int16(cor2_sar_isp_surf);
hn_mean_ku = uint32(mean_cai_fai_sar_surf);
pri_lrm_l1b_ku = uint32(pri_sar_isp_surf .* T0_sar_surf_nadir * pri_T0_unit_conv_chd * 1e12);

%----------M. Waveform related variables ----------------------------------
waveform_scale_factor_ku = zeros(N_total_surf_loc,1);
sar_power_waveform_ku = zeros(N_total_surf_loc, N_samples_sar_chd * zp_fact_range_cnf);

if(no_amb==0)
    for i_surf=1:N_total_surf_loc
        waveform_scale_factor_ku(i_surf) = single(max(wfm_cor_i2q2_sar_ku(i_surf,:))*10^0.3 / (2^16-1));
        sar_power_waveform_ku(i_surf,:) = uint16(round(wfm_cor_i2q2_sar_ku(i_surf,:).*10^0.3 ./ waveform_scale_factor_ku(i_surf)));
    end
elseif(no_amb==1)
   for i_surf=1:N_total_surf_loc
        waveform_scale_factor_ku(i_surf) = single(max(wfm_cor_i2q2_sar_ku_noDopp(i_surf,:))*10^0.3 / (2^16-1));
        sar_power_waveform_ku(i_surf,:) = uint16(round(wfm_cor_i2q2_sar_ku_noDopp(i_surf,:).*10^0.3 ./ waveform_scale_factor_ku(i_surf)));
    end     
end
snr_estimation_ku = zeros(N_total_surf_loc,1); % TBD
sigma0_scaling_factor_ku = int16(wfm_scaling_factor_sar_ku * 1e2);

%----------P. Stack characterisation variables	-----------------------------
look_angle_start_ku = zeros(N_total_surf_loc,1);
look_angle_stop_ku = zeros(N_total_surf_loc,1);
doppler_angle_start_ku = zeros(N_total_surf_loc,1);
doppler_angle_stop_ku = zeros(N_total_surf_loc,1);
pointing_angle_start_ku = zeros(N_total_surf_loc,1);
pointing_angle_stop_ku = zeros(N_total_surf_loc,1);


%finding indexes and number of contributing beams
% for i_surf=1:N_total_surf_loc
%     N_beams_contributing_noamb(i_surf) = max(N_beams_stack)-length(find(stack_mask_vector(i_surf,:)==1));
%     N_beams_contributing_amb(i_surf) = max(N_beams_stack)-length(find(geocorr_mask_vector(i_surf,:)==1));
%     first_contrib_beam(i_surf)=(find(geocorr_mask_vector(i_surf,:)==1));
% end


num_total_looks_ku = uint16(N_beams_stack');
num_contributing_looks_ku= num_total_looks_ku;
% 
% 
% if(no_amb==0)
%    num_contributing_looks_ku = uint16(N_beams_contributing_amb'); 
% elseif(no_amb==1)
%    num_contributing_looks_ku = uint16(N_beams_contributing_noamb'); 
% end
% first_beam_index = num_total_looks_ku-num_contributing_looks_ku;

for i_surf=1:N_total_surf_loc
    look_angle_start_ku(i_surf) = int16(look_ang_surf(i_surf,1)'*1e6);
    look_angle_stop_ku(i_surf) = int16(look_ang_surf(i_surf,N_beams_stack(i_surf))*1e6);
    doppler_angle_start_ku(i_surf) = int16((doppler_ang_surf(i_surf,1))'*1e6);
    doppler_angle_stop_ku (i_surf) = int16((doppler_ang_surf(i_surf,N_beams_stack(i_surf)))'*1e6);
    pointing_angle_start_ku(i_surf) = int16(pointing_ang_surf(i_surf,1)*1e6)';
    pointing_angle_stop_ku (i_surf) = int16(pointing_ang_surf(i_surf,N_beams_stack(i_surf))'*1e6);
end


gaussian_fitting_skewness_ku = int16(stack_skewness' * 1e2);
gaussian_fitting_kurtosis_ku = int16(stack_kurtosis' * 1e2);
gaussian_fitting_std_ku = uint16(stack_std' * 1e6);
gaussian_fitting_centre_look_ku = int16(stack_look_ang_centre' * 1e6);
gaussian_fitting_centre_pointing_ku = int16(stack_pointing_ang_centre' * 1e6);

%% PACKING L1B

%----------A. Time variables ----------------------------------------------
l1_mode_id_ku_name = 'l1_mode_id_ku';
nccreate(filename,l1_mode_id_ku_name,format_key,netcdf_v4_format,data_type_key,uint8_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite( filename,l1_mode_id_ku_name,l1_mode_id_ku);
ncwriteatt(filename,['/' l1_mode_id_ku_name],long_name_att,'L1 mode ID (Ku-band)');
ncwriteatt(filename,['/' l1_mode_id_ku_name],comment_att,'L1 Mode Identifier. Each L1 record type has a unique ID, as follows: 0 = Null record, 1 = LRM L1b Ku, 2 = LRM L1b C, 3 = SAR Raw L1A Ku, 4 = SAR Raw L1B-S Ku, 5 = SAR Raw L1B Ku, 6 = SAR RMC L1A Ku, 7 = SAR RMC L1B-S Ku, 8 = SAR RMC L1B Ku, 9 = CAL1 LRM L1b Ku, 10 = CAL1 LRM L1b C, 11 = CAL1 SAR L1b Ku, 12 = CAL1 SAR L1b C, 13 = CAL1 INSTR L1b Ku, 14 = CAL1 INSTR L1b C, 15 = Autocal L1b, 16 = CAL2 L1b Ku, 17 = CAL2 L1b C');

time_day_ku_name = 'time_day_ku';
nccreate(filename,time_day_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite( filename,time_day_ku_name,time_day_ku);
ncwriteatt(filename,['/' time_day_ku_name],long_name_att,'Day from 1 Jan 2000 (Ku-band)');
ncwriteatt(filename,['/' time_day_ku_name],units_att,day_units);
ncwriteatt(filename,['/' time_day_ku_name],comment_att,'It contains the day from 1 Jan 2000. Time is in TAI. Time refers to the instant the L1B waveform (which is a combination of pulses for both LRM and SAR) touches the surface');

time_seconds_ku_name = 'time_seconds_ku';
nccreate(filename,time_seconds_ku_name,format_key,netcdf_v4_format,data_type_key,uint64_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite( filename,time_seconds_ku_name,time_seconds_ku);
ncwriteatt(filename,['/' time_seconds_ku_name],long_name_att,'Seconds of the day, with microsecond resolution (Ku-band)');
ncwriteatt(filename,['/' time_seconds_ku_name],units_att,seconds_units);
ncwriteatt(filename,['/' time_seconds_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' time_seconds_ku_name],comment_att,'Second of the day, with a scaling factor of 1.0e-6, providing a microsecond resolution. Time is in TAI. Time refers to the instant the L1B waveform (which is a combination of pulses for both LRM and SAR) touches the surface');

tm_source_sequence_counter_ku_name = 'tm_source_sequence_counter_ku';
nccreate(filename,tm_source_sequence_counter_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,tm_source_sequence_counter_ku_name,tm_source_sequence_counter_ku);
ncwriteatt(filename,['/' tm_source_sequence_counter_ku_name],long_name_att,'Instrument source sequence counter (Ku-band)');
ncwriteatt(filename,['/' tm_source_sequence_counter_ku_name],units_att,number_units);
ncwriteatt(filename,['/' tm_source_sequence_counter_ku_name],comment_att,'This is the instrument record counter, copied from ISP');

l1b_record_counter_ku_name = 'l1b_record_counter_ku';
nccreate(filename,l1b_record_counter_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,l1b_record_counter_ku_name,l1b_record_counter_ku);
ncwriteatt(filename,['/' l1b_record_counter_ku_name],long_name_att,'L1B record counter (Ku-band)');
ncwriteatt(filename,['/' l1b_record_counter_ku_name],units_att,number_units);
ncwriteatt(filename,['/' l1b_record_counter_ku_name],comment_att,' The L1B record counter, starting from 0 for each product');

%----------B. Orbit and attitude variables ------------------------------
latitude_ku_name = 'latitude_ku';
nccreate(filename,latitude_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,latitude_ku_name,latitude_ku);
ncwriteatt(filename,['/' latitude_ku_name],long_name_att,'latitude (positive N, negative S) (Ku-band)');
ncwriteatt(filename,['/' latitude_ku_name],units_att,degrees_units);
ncwriteatt(filename,['/' latitude_ku_name],scale_factor_att,'1e-7');
ncwriteatt(filename,['/' latitude_ku_name],comment_att,'Latitude of measurement [-90, +90]: Positive at Nord, Negative at South');

longitude_ku_name = 'longitude_ku';
nccreate(filename,longitude_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,longitude_ku_name,longitude_ku);
ncwriteatt(filename,['/' longitude_ku_name],long_name_att,'Longitude (degrees East) (Ku-band)');
ncwriteatt(filename,['/' longitude_ku_name],units_att,degrees_units);
ncwriteatt(filename,['/' longitude_ku_name],scale_factor_att,'1e-7');
ncwriteatt(filename,['/' longitude_ku_name],comment_att,'Longitude of measurement [0, 360]: Positive at East');

com_altitude_ku_name = 'com_altitude_ku';
nccreate(filename,com_altitude_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,com_altitude_ku_name,com_altitude_ku);
ncwriteatt(filename,['/' com_altitude_ku_name],long_name_att,'CoM altitude (Ku-band)');
ncwriteatt(filename,['/' com_altitude_ku_name],units_att,meters_units);
ncwriteatt(filename,['/' com_altitude_ku_name],add_offset_att,'1.3e6');
ncwriteatt(filename,['/' com_altitude_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' com_altitude_ku_name],comment_att,'Altitude of the satellite Centre of Mass');

com_altitude_rate_ku_name = 'com_altitude_rate_ku';
nccreate(filename,com_altitude_rate_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,com_altitude_rate_ku_name,com_altitude_rate_ku);
ncwriteatt(filename,['/' com_altitude_rate_ku_name],long_name_att,'CoM altitude rate (Ku-band)');
ncwriteatt(filename,['/' com_altitude_rate_ku_name],units_att,meters_per_second_units);
ncwriteatt(filename,['/' com_altitude_rate_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' com_altitude_rate_ku_name],comment_att,'Instantaneous altitude rate at the Centre of Mass');

com_velocity_vector_ku_name = 'com_velocity_vector_ku';
nccreate(filename,com_velocity_vector_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{space_3D_dimension ,space_3D_dimension_size,ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,com_velocity_vector_ku_name,com_velocity_vector_ku);
ncwriteatt(filename,['/' com_velocity_vector_ku_name],long_name_att,'Velocity vector in ITRF, components: [1] x, [2] y, [3] z (Ku-band)');
ncwriteatt(filename,['/' com_velocity_vector_ku_name],units_att,meters_per_second_units);
ncwriteatt(filename,['/' com_velocity_vector_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' com_velocity_vector_ku_name],comment_att,'Velocity Vector at the centre of mass in ITRF. The 3 components are given according to the ''space_3D'' dimension: [1] x, [2] y, [3] z');

satellite_mispointing_ku_name = 'satellite_mispointing_ku';
nccreate(filename,satellite_mispointing_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{space_3D_dimension ,space_3D_dimension_size,ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,satellite_mispointing_ku_name,satellite_mispointing_ku);
ncwriteatt(filename,['/' satellite_mispointing_ku_name],long_name_att,'Mispointing angle, measures by STRs: [1] Roll, [2] Pitch, [3] Yaw (Ku-band)');
ncwriteatt(filename,['/' satellite_mispointing_ku_name],units_att,degrees_units);
ncwriteatt(filename,['/' satellite_mispointing_ku_name],scale_factor_att,'1e-7');
ncwriteatt(filename,['/' satellite_mispointing_ku_name],comment_att,'Attitude mispointing, measured by STRs and post-processed by AOCS or by ground facility. The 3 components are given according to the ''space_3D'' dimension: [1] Roll, [2] Pitch, [3] Yaw. This variable includes the "mispointing bias" given by the variable mispointing_bias_ku. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: roll = pitch = yaw = 0');

mispointing_bias_ku_name = 'mispointing_bias_ku';
nccreate(filename,mispointing_bias_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{space_3D_dimension, space_3D_dimension_size});
ncwrite(filename,mispointing_bias_ku_name,mispointing_bias_ku);
ncwriteatt(filename,['/' mispointing_bias_ku_name],long_name_att,'Mispointing bias from in-flight calibration: [1] Roll, [2] Pitch, [3] Yaw (Ku-band)');
ncwriteatt(filename,['/' mispointing_bias_ku_name],units_att,degrees_units);
ncwriteatt(filename,['/' mispointing_bias_ku_name],scale_factor_att,'1e-7');
ncwriteatt(filename,['/' mispointing_bias_ku_name],comment_att,'The mispointing bias will correct for structural misalignment not accounted for during on-ground characterization. It will provide absolute mispointing value. It will be estimated using altimeter''s measurements data. The mispointing bias is given for [1] roll, [2] pitch and [3] yaw and it is constant within the product.');

%----------C. Configuration and quality variables -------------------------
l1_instrument_configuration_ku_name = 'l1_instrument_configuration_ku';
nccreate(filename,l1_instrument_configuration_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,l1_instrument_configuration_ku_name,l1_instrument_configuration_ku_surf);
ncwriteatt(filename,['/' l1_instrument_configuration_ku_name],long_name_att,'Instrument configuration L1 (Ku-band)');
ncwriteatt(filename,['/' l1_instrument_configuration_ku_name],comment_att,'It contains flags from ISPs related to P4 configuration and tracking. For more details, see ESA IODD "JC-ID-ESA-GP-0059"');
% 
l1b_mcd_ku_name = 'l1b_mcd_ku';
nccreate(filename,l1b_mcd_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
% ncwrite(filename,l1b_mcd_ku_name,l1b_mcd_ku);
ncwriteatt(filename,['/' l1b_mcd_ku_name],long_name_att,'Measurement confidence data (MCD) L1A (Kuband)');
ncwriteatt(filename,['/' l1b_mcd_ku_name],comment_att,'It contains flags related to the L1A processing chain. For more details, see ESA IODD "JC-ID-ESA-GP-0059"');

%----------D. Altimeter range variables -----------------------------------  
altimeter_range_calibrated_ku_name = 'altimeter_range_calibrated_ku';
nccreate(filename,altimeter_range_calibrated_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,altimeter_range_calibrated_ku_name,altimeter_range_calibrated_ku);
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],long_name_att,'Calibrated 1-way range: CoM to middle range window (at sample Ns/2 from 0) (Ku-band)');
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],units_att,meters_units);
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],add_offset_att,'1.3e6');
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],comment_att,'This is the 1-way distance from the satellite''s Center of Mass to the middle of the range window (sample Ns/2 from 0). It includes the following range calibrations: (a) range_corr_internal_delay, (b) range_corr_com, (c) range_corr_doppler [LRM only]. Note: the actual altimeter clock (variable ''altimeter_clock_ku'') has been used to compute the altimeter range');

range_corr_internal_delay_ku_name = 'range_corr_internal_delay_ku';
nccreate(filename,range_corr_internal_delay_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,range_corr_internal_delay_ku_name,range_corr_internal_delay_ku);
ncwriteatt(filename,['/' range_corr_internal_delay_ku_name],long_name_att,'Range correction:internal path delay (waveguides + CAL1) (Ku-band)');
ncwriteatt(filename,['/' range_corr_internal_delay_ku_name],units_att,meters_units);
ncwriteatt(filename,['/' range_corr_internal_delay_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' range_corr_internal_delay_ku_name],comment_att,'The internal path delay is the sum of 2 components: (a) waveguides component: constant part measured pre-launch, (b) CAL1 P4 internal delay: periodically measured on flight');

range_corr_com_ku_name = 'range_corr_com_ku';
nccreate(filename,range_corr_com_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,range_corr_com_ku_name,range_corr_com_ku);
ncwriteatt(filename,['/' range_corr_com_ku_name],long_name_att,'Range correction: distance CoM antenna(Ku-band)');
ncwriteatt(filename,['/' range_corr_com_ku_name],units_att,meters_units);
ncwriteatt(filename,['/' range_corr_com_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' range_corr_com_ku_name],comment_att,'It is the distance in the z-component between the Centre of Mass and the Antenna reference point. The Centre of Mass is assumed constant in the GPP, although it is likely to vary in P4');

%----------E. Altimeter power variables -----------------------------------
attenuator_calibrated_ku_name = 'attenuator_calibrated_ku';
nccreate(filename,attenuator_calibrated_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,attenuator_calibrated_ku_name,attenuator_calibrated_ku);
ncwriteatt(filename,['/' attenuator_calibrated_ku_name],long_name_att,'Power scaling: ATT calibrated (from ATT table) (Ku-band)');
ncwriteatt(filename,['/' attenuator_calibrated_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' attenuator_calibrated_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' attenuator_calibrated_ku_name],comment_att,'The ISPs contain the ATT command applied within a radar cycle. The actual ATT value, that slightly differs from the nominal one, is extracted from a look-up table');

altimeter_power_drift_ku_name = 'altimeter_power_drift_ku';
nccreate(filename,altimeter_power_drift_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,altimeter_power_drift_ku_name,altimeter_power_drift_ku);
ncwriteatt(filename,['/' altimeter_power_drift_ku_name],long_name_att,'Power scaling: power instrument aging correction (CAL1) (Ku-band)');
ncwriteatt(filename,['/' altimeter_power_drift_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' altimeter_power_drift_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' altimeter_power_drift_ku_name],comment_att,'This is a measure of the instrument power decay over the time, regularly measured from CAL1 impulse response with respect to a reference at the beginning of the mission');

power_corr_digital_processing_ku_name = 'power_corr_digital_processing_ku';
nccreate(filename,power_corr_digital_processing_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type);
ncwrite(filename,power_corr_digital_processing_ku_name,power_corr_digital_processing_ku);
ncwriteatt(filename,['/' power_corr_digital_processing_ku_name],long_name_att,'Power scaling: onboard digital processing (Ku-band)');
ncwriteatt(filename,['/' power_corr_digital_processing_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' power_corr_digital_processing_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' power_corr_digital_processing_ku_name],comment_att,'Each measurement mode (LRM, SAR, RMC) is subject to a specific on-board digital processing that affects the signal scaling. This factor includes all the digital processing scaling factors, which are constant for each mode');

power_scaling_to_antenna_ku_name = 'power_scaling_to_antenna_ku';
nccreate(filename,power_scaling_to_antenna_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,power_scaling_to_antenna_ku_name,power_scaling_to_antenna_ku);
ncwriteatt(filename,['/' power_scaling_to_antenna_ku_name],long_name_att,'Power scaling: overall gain from raw waveform to antenna flange (Ku-band)');
ncwriteatt(filename,['/' power_scaling_to_antenna_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' power_scaling_to_antenna_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' power_scaling_to_antenna_ku_name],comment_att,'The overall instrument gain, including: (a) attenuator_calibrated, (b) altimeter_power_drift, (c) power_corr_digital_processing, (d) RFU gains. Antenna gain is excluded. This variable is applied at L1B for LRM and at L1A for SAR to calibrate the power level. It follows that the L1B waveform power is the received power at the antenna flange');

%----------F. Altimeter engineering variables -----------------------------
altimeter_clock_ku_name = 'altimeter_clock_ku';
nccreate(filename,altimeter_clock_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,altimeter_clock_ku_name,altimeter_clock_ku);
ncwriteatt(filename,['/' altimeter_clock_ku_name],long_name_att,'Altimeter clock (Ku-band)');
ncwriteatt(filename,['/' altimeter_clock_ku_name],units_att,Hz_units);
ncwriteatt(filename,['/' altimeter_clock_ku_name],add_offset_att,'395e6');
ncwriteatt(filename,['/' altimeter_clock_ku_name],scale_factor_att,'1e-9');
ncwriteatt(filename,['/' altimeter_clock_ku_name],comment_att,'This is the actual altimeter clock. The altimeter clock is based upon the USO clock. The nominal USO clock is 10MHz, while the nominal altimeter clock is 395 MHz. The actual USO clock is provided regularly as a drift information: based on the USO drift, the actual altimeter clock is computed, assuming a linear dependency');

tm_h0_ku_name = 'tm_h0_ku';
nccreate(filename,tm_h0_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,tm_h0_ku_name,tm_h0_ku);
ncwriteatt(filename,['/' tm_h0_ku_name],long_name_att,'H0 telemetry (same in 1 radar cycle) (Ku-band)');
ncwriteatt(filename,['/' tm_h0_ku_name],units_att,T0d64_units);
ncwriteatt(filename,['/' tm_h0_ku_name],comment_att,'The H0 (initial altitude instruction) value, copied from ISP. The altimeter provides one H0 value per radar cycle. It is used in combination with the COR2 to compute the CAI and FAI instructions. NOTE: T0 = 1/altimeter_clock_ku');

tm_cor2_ku_name = 'tm_cor2_ku';
nccreate(filename,tm_cor2_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,tm_cor2_ku_name,tm_cor2_ku);
ncwriteatt(filename,['/' tm_cor2_ku_name],long_name_att,'COR2 telemetry (same in 1 radar cycle) (Ku-band)');
ncwriteatt(filename,['/' tm_cor2_ku_name],units_att,T0d16d64_units);
ncwriteatt(filename,['/' tm_cor2_ku_name],comment_att,'The COR2 (altitude rate estimation) value, copied from ISP. The altimeter estimates one COR2 value per radar cycle. It is used in combination with the H0 to compute the CAI and FAI instructions. NOTE: T0 = 1/altimeter_clock_ku');

hn_mean_ku_name = 'hn_mean_ku';
nccreate(filename,hn_mean_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,hn_mean_ku_name,hn_mean_ku);
ncwriteatt(filename,['/' hn_mean_ku_name],long_name_att,'Mean altitude instruction,	for	monitoring purpose (Ku-band)');
ncwriteatt(filename,['/' hn_mean_ku_name],units_att,T0d16d64_units);
ncwriteatt(filename,['/' hn_mean_ku_name],comment_att,'This field is for monitoring the H0/COR2 decoding into CAI (Coarse Altitude instruction) and FAI (Fine Altitude Instruction). In LRM it is the altitude instruction (CAI+FAI) averaged and rounded over the radar cycle. In SAR, it is the altitude instruction (CAI+FAI) averaged and rounded over the burst closest to the nadir of the surface location. NOTE: T0 = 1/altimeter_clock_ku');

pulse_repetition_interval_ku_name = 'pulse_repetition_interval_ku';
nccreate(filename,pulse_repetition_interval_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,pulse_repetition_interval_ku_name,pri_lrm_l1b_ku);
ncwriteatt(filename,['/' pulse_repetition_interval_ku_name],long_name_att,'PRI converted into seconds (Ku-band)');
ncwriteatt(filename,['/' pulse_repetition_interval_ku_name],units_att,seconds_units);
ncwriteatt(filename,['/' pulse_repetition_interval_ku_name],scale_factor_att,'1e-12');
ncwriteatt(filename,['/' pulse_repetition_interval_ku_name],comment_att,'The ''Pulse Repetition Interval''. PRI is constant within all received pulses in a radar cycle, but it can change within consecutive radar cycles. It is provided in counters of (T0*8) by the altimeter and converted in seconds by the L1 processor, with T0 = 1/altimeter_clock_ku');

%----------M. Waveform related variables -----------------------------
sar_power_waveform_ku_name = 'sar_power_waveform_ku';
nccreate(filename,sar_power_waveform_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ns_dimension, N_samples_sar_chd*zp_fact_range_cnf,ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,sar_power_waveform_ku_name,sar_power_waveform_ku.');
ncwriteatt(filename,['/' sar_power_waveform_ku_name],long_name_att,'SAR	Power Echo waveform: scaled	0-65535 (Ku-band)');
ncwriteatt(filename,['/' sar_power_waveform_ku_name],units_att,number_units);
ncwriteatt(filename,['/' sar_power_waveform_ku_name],comment_att,'The SAR L1B Power waveforms is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B_MCD flags, (c) final scaling, given in the variable ''waveform_scale_factor'', in order to best fit the waveform into 2 bytes');

stack_mask_range_bin_ku_name = 'stack_mask_range_bin_ku';
nccreate(filename,stack_mask_range_bin_ku_name,format_key,netcdf_v4_format,data_type_key,uint8_type,dimensions_key,{nl_dimension, nl_dimension_size,ku_rec_dimension, N_total_surf_loc});

if(no_amb==1)
    ncwrite(filename,stack_mask_range_bin_ku_name,uint8(-1+ceil(stack_mask_vector/zp_fact_range_cnf)).');
elseif(no_amb==0)
    ncwrite(filename,stack_mask_range_bin_ku_name,uint8(-1+ceil(geocorr_mask_vector/zp_fact_range_cnf)).');
end
ncwriteatt(filename,['/' stack_mask_range_bin_ku_name],long_name_att,'Range bin stack mask ');
ncwriteatt(filename,['/' stack_mask_range_bin_ku_name],units_att,number_units);
ncwriteatt(filename,['/' stack_mask_range_bin_ku_name],scale_factor_att,zp_fact_range_cnf);
ncwriteatt(filename,['/' stack_mask_range_bin_ku_name],comment_att,'Before the stack is multi-looked, the different looks are cropped according to these values, in order to remove the ambiguities. For each look, the number of the first cropped sample is provided. The scale factor is equal to the range_oversampling_factor');


waveform_scale_factor_ku_name = 'waveform_scale_factor_ku';
nccreate(filename,waveform_scale_factor_ku_name,format_key,netcdf_v4_format,data_type_key,float_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,waveform_scale_factor_ku_name,waveform_scale_factor_ku);
ncwriteatt(filename,['/' waveform_scale_factor_ku_name],long_name_att,'Echo Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
ncwriteatt(filename,['/' waveform_scale_factor_ku_name],units_att,W_per_count_units);
ncwriteatt(filename,['/' waveform_scale_factor_ku_name],comment_att,'The L1B waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = power_waveform(ku_rec, Ns) * waveform_scale_factor(ku_rec)');

zero_padding_ku_name = 'range_oversampling_factor_ku';
nccreate(filename,zero_padding_ku_name,format_key,netcdf_v4_format,data_type_key,int8_type);
ncwrite(filename,zero_padding_ku_name,zp_fact_range_cnf);
ncwriteatt(filename,['/' zero_padding_ku_name],long_name_att,'Oversampling factor used in the range compression (FFT)');
ncwriteatt(filename,['/' zero_padding_ku_name],units_att,number_units);
ncwriteatt(filename,['/' zero_padding_ku_name],comment_att,'The instrument samples the waveforms with a 395 MHz clock, providing a nominal_sampling = c / 395e6 / 2 = ~0.379m (with c=speed of light). In addition, the ground processor can apply an oversampling factor, providing a waveform_sampling = nominal_sampling / range_oversampling_factor. Note that the altimeter range resolution is fixed and given by the chirp bandwidth of 320 MHz: c / 320e6 / 2 = ~0.468m');

% 
snr_estimation_ku_name = 'snr_estimation_ku';
nccreate(filename,snr_estimation_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
% ncwrite(filename,snr_estimation_ku_name,snr_estimation_ku);
ncwriteatt(filename,['/' snr_estimation_ku_name],long_name_att,'SNR estimation on	L1B	waveform');
ncwriteatt(filename,['/' snr_estimation_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' snr_estimation_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' snr_estimation_ku_name],comment_att,'SNR estimation done on the L1B waveform as: peak / noise_level. Where noise_level is estimated from the first n (TBC) samples');


sigma0_scaling_factor_ku_name = 'sigma0_scaling_factor_ku';
nccreate(filename,sigma0_scaling_factor_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,sigma0_scaling_factor_ku_name,sigma0_scaling_factor_ku);
ncwriteatt(filename,['/' sigma0_scaling_factor_ku_name],long_name_att,'Power scaling: Scaling factor for sigma-0 evaluation (from antenna flange)');
ncwriteatt(filename,['/' sigma0_scaling_factor_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' sigma0_scaling_factor_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' sigma0_scaling_factor_ku_name],comment_att,'This is a scaling factor in order to retrieve sigma-0 from the L1B waveform. It includes antenna gains and geometry satellite - surface. It is not applied to the L1B waveforms');

%----------P. Stack characterisation variables	-----------------------------
num_total_looks_ku_name = 'num_total_looks_ku';
nccreate(filename,num_total_looks_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,num_total_looks_ku_name,num_total_looks_ku);
ncwriteatt(filename,['/' num_total_looks_ku_name],long_name_att,'Number of total looks in the stack (monitoring purpose) (Ku-band)');
ncwriteatt(filename,['/' num_total_looks_ku_name],units_att,number_units);
ncwriteatt(filename,['/' num_total_looks_ku_name],comment_att,'The total number of looks available within the stack. Some of these looks may be weighted or removed prior to multi-looking');

num_contributing_looks_ku_name = 'num_contributing_looks_ku';
nccreate(filename,num_contributing_looks_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,num_contributing_looks_ku_name,num_contributing_looks_ku);
ncwriteatt(filename,['/' num_contributing_looks_ku_name],long_name_att,'Number of contributing looks in the L1B wavefrom (non-0 weight) (Ku-band)');
ncwriteatt(filename,['/' num_contributing_looks_ku_name],units_att,number_units);
ncwriteatt(filename,['/' num_contributing_looks_ku_name],comment_att,'The number of looks that have been actually used within the multi-looking. The number of looks with a non-0 weight');

look_angle_start_ku_name = 'look_angle_start_ku';
nccreate(filename,look_angle_start_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,look_angle_start_ku_name,look_angle_start_ku);
ncwriteatt(filename,['/' look_angle_start_ku_name],long_name_att,'Angle of first look  (Ku-band)');
ncwriteatt(filename,['/' look_angle_start_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' look_angle_start_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' look_angle_start_ku_name],comment_att,'Look angle of the first contributing look (non-0 weight) to the L1B waveform');

look_angle_stop_ku_name = 'look_angle_stop_ku';
nccreate(filename,look_angle_stop_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,look_angle_stop_ku_name,look_angle_stop_ku);
ncwriteatt(filename,['/' look_angle_stop_ku_name],long_name_att,'Angle of last look  (Ku-band)');
ncwriteatt(filename,['/' look_angle_stop_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' look_angle_stop_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' look_angle_stop_ku_name],comment_att,'Look angle of the last contributing look (non-0 weight) to the L1B waveform');

doppler_angle_start_ku_name = 'doppler_angle_start_ku';
nccreate(filename,doppler_angle_start_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,doppler_angle_start_ku_name,doppler_angle_start_ku);
ncwriteatt(filename,['/' doppler_angle_start_ku_name],long_name_att,'Angle of first look  (Ku-band)');
ncwriteatt(filename,['/' doppler_angle_start_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' doppler_angle_start_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' doppler_angle_start_ku_name],comment_att,'Doppler angle of the first contributing look (non-0 weight) to the L1B waveform');
 
doppler_angle_stop_ku_name = 'doppler_angle_stop_ku';
nccreate(filename,doppler_angle_stop_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,doppler_angle_stop_ku_name,doppler_angle_stop_ku);
ncwriteatt(filename,['/' doppler_angle_stop_ku_name],long_name_att,'Angle of last contributing look  (Ku-band)');
ncwriteatt(filename,['/' doppler_angle_stop_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' doppler_angle_stop_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' doppler_angle_stop_ku_name],comment_att,'Doppler angle of the last contributing look to the L1B waveform');

pointing_angle_start_ku_name = 'pointing_angle_start_ku';
nccreate(filename,pointing_angle_start_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,pointing_angle_start_ku_name,pointing_angle_start_ku);
ncwriteatt(filename,['/' pointing_angle_start_ku_name],long_name_att,'Angle of first contributing look  (Ku-band)');
ncwriteatt(filename,['/' pointing_angle_start_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' pointing_angle_start_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' pointing_angle_start_ku_name],comment_att,'Pointing angle of the first contributing look to the L1B waveform');
 
pointing_angle_stop_ku_name = 'pointing_angle_stop_ku';
nccreate(filename,pointing_angle_stop_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,pointing_angle_stop_ku_name,pointing_angle_stop_ku);
ncwriteatt(filename,['/' pointing_angle_stop_ku_name],long_name_att,'Angle of last contributing look  (Ku-band)');
ncwriteatt(filename,['/' pointing_angle_stop_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' pointing_angle_stop_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' pointing_angle_stop_ku_name],comment_att,'Pointing angle of the last contributing look to the L1B waveform');

gaussian_fitting_skewness_ku_name = 'gaussian_fitting_skewness_ku';
nccreate(filename,gaussian_fitting_skewness_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,gaussian_fitting_skewness_ku_name,gaussian_fitting_skewness_ku);
ncwriteatt(filename,['/' gaussian_fitting_skewness_ku_name],long_name_att,'Gaussian Power fitting: skewness (Ku-band)');
ncwriteatt(filename,['/' gaussian_fitting_skewness_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' gaussian_fitting_skewness_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' gaussian_fitting_skewness_ku_name],comment_att,'Skewness of the Gaussian that fits the integrated power of the looks within a stack. The skewness indicates how symmetric or asymmetric the power within the stack is');

gaussian_fitting_kurtosis_ku_name = 'gaussian_fitting_kurtosis_ku';
nccreate(filename,gaussian_fitting_kurtosis_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,gaussian_fitting_kurtosis_ku_name,gaussian_fitting_kurtosis_ku);
ncwriteatt(filename,['/' gaussian_fitting_kurtosis_ku_name],long_name_att,'Gaussian Power fitting: kurtosis (Ku-band)');
ncwriteatt(filename,['/' gaussian_fitting_kurtosis_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' gaussian_fitting_kurtosis_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' gaussian_fitting_kurtosis_ku_name],comment_att,'Kurtosis of the Gaussian that fits the integrated power of the looks within a stack. Kurtosis is a measure of peakiness');

gaussian_fitting_std_ku_name = 'gaussian_fitting_std_ku';
nccreate(filename,gaussian_fitting_std_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,gaussian_fitting_std_ku_name,gaussian_fitting_std_ku);
ncwriteatt(filename,['/' gaussian_fitting_std_ku_name],long_name_att,'Gaussian Power fitting: STD wrt look angle (Ku-band)');
ncwriteatt(filename,['/' gaussian_fitting_std_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' gaussian_fitting_std_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' gaussian_fitting_std_ku_name],comment_att,'Standard deviation of the Gaussian that fits the integrated power of the looks within a stack. It is given with respect to the look angle. The width at -3dB of this Gaussian can be retrieved the following way: width_3db = 2*sqrt(2*ln2)*gaussian_fitting_std');

gaussian_fitting_centre_look_ku_name = 'gaussian_fitting_centre_look_ku';
nccreate(filename,gaussian_fitting_centre_look_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,gaussian_fitting_centre_look_ku_name,gaussian_fitting_centre_look_ku);
ncwriteatt(filename,['/' gaussian_fitting_centre_look_ku_name],long_name_att,'Gaussian Power fitting: centre wrt look angle (Ku-band)');
ncwriteatt(filename,['/' gaussian_fitting_centre_look_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' gaussian_fitting_centre_look_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' gaussian_fitting_centre_look_ku_name],comment_att,'Position of the center of the Gaussian that fits the integrated power of the looks within a stack, with respect to look angle');

gaussian_fitting_centre_pointing_ku_name = 'gaussian_fitting_centre_pointing_ku';
nccreate(filename,gaussian_fitting_centre_pointing_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_surf_loc});
ncwrite(filename,gaussian_fitting_centre_pointing_ku_name,gaussian_fitting_centre_pointing_ku);
ncwriteatt(filename,['/' gaussian_fitting_centre_pointing_ku_name],long_name_att,'Gaussian Power fitting: centre wrt pointing angle (Ku-band)');
ncwriteatt(filename,['/' gaussian_fitting_centre_pointing_ku_name],units_att,rad_units);
ncwriteatt(filename,['/' gaussian_fitting_centre_pointing_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' gaussian_fitting_centre_pointing_ku_name],comment_att,'Position of the centre of the Gaussian that fits the integrated power of the looks within a stack, with respect to pointing angle');

ncwriteatt(filename,'/','creation_time',datestr(now));
ncwriteatt(filename,'/','data_info','data simulated by ESA and processed by isardSAT');

time = toc(t6);
minutes_reading = floor(time/60);
secs_reading = time - minutes_reading*60;
disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
