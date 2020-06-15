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
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/07/2014)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils


t6 = tic;

if strcmp(TEST,'S6_C1A_SAR')
    filename = 'S6_P4_SIM_RAW_L1A__20190119T064000_20190119T064019_T01.nc';
elseif strcmp(TEST,'S6_C1A_RMC')
    filename = 'S6_P4_SIM_RMC_L1A__20190119T064000_20190119T064019_T01.nc';    
elseif strcmp(TEST,'S6_C1B_SAR')
    filename = 'S6_P4_SIM_RAW_L1A__20210929T064000_20210929T064019_T02.nc';    
elseif strcmp(TEST,'S6_C1B_RMC')
    filename = 'S6_P4_SIM_RMC_L1A__20210929T064000_20210929T064019_T02.nc';    
elseif strcmp(TEST,'S6_C1C_SAR')
    filename = 'S6_P4_SIM_RAW_L1A__20210929T064000_20210929T064019_T03.nc';    
elseif strcmp(TEST,'S6_C1C_RMC')
    filename = 'S6_P4_SIM_RMC_L1A__20210929T064000_20210929T064019_T03.nc';    
end



dimensions_key = 'Dimensions';
format_key = 'Format';
data_type_key = 'DataType';

netcdf_v4_format = 'netcdf4';

ku_rec_dimension = 'Ku_rec';
np_dimension = 'Np';
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
T0x4_units = 'T0*4';
T0d64_units = 'T0/64';
T0x8_units = 'T0*8';
T0d16d64_units = 'T0/16/64';
sqrtW_per_count_units = 'sqrt(Watt)/#';

%% CODING L1A
%----------A. Time variables ----------------------------------------------
if (Process_ID2 == 58) % TM_ECHO_SAR
    l1_mode_id_ku = uint8(zeros(N_total_bursts_sar_ku_isp,1)+3); % SAR Raw L1A Ku
elseif (Process_ID2 == 59) % TM_ECHO_RMC
    l1_mode_id_ku = uint8(zeros(N_total_bursts_sar_ku_isp,1)+6); % SAR RMC L1A Ku
else % unknown
    l1_mode_id_ku = uint8(zeros(N_total_bursts_sar_ku_isp,1)+0);
end
time_day_ku = uint32(time_sar_ku / sec_in_day_cst);
time_seconds_ku = uint64((time_sar_ku - double(time_day_ku) * sec_in_day_cst) * 1e6);
source_seq_count_sar_isp= 0:(N_total_bursts_sar_ku_isp-1);
tm_source_sequence_counter_ku = uint16(source_seq_count_sar_isp).';
burst_counter_ku = uint16(rec_counter_sar_isp-1);

%----------B. Orbit and attitude variables ------------------------------
latitude_ku = int32(lat_sar_sat * 1e7);
for i_lon = 1:size(lon_sar_sat,1)
    if (lon_sar_sat(i_lon) > 180)
        longitude_ku = uint32((lon_sar_sat +360)* 1e7); 
    else
        longitude_ku = uint32((lon_sar_sat )* 1e7); 
    end
end
com_altitude_ku = uint32((alt_sar_sat - 1.3e6) * 1e4);
com_altitude_rate_ku = int32(alt_rate_sar_sat * 1e4);
com_velocity_vector_ku = int32([x_vel_sat_sar * 1e4; y_vel_sat_sar * 1e4; z_vel_sat_sar * 1e4]);
satellite_mispointing_ku = int32([pitch_sar * 180/pi_cst * 1e7; roll_sar * 180/pi_cst * 1e7; yaw_sar * 180/pi_cst * 1e7]);
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
%     l1_instrument_configuration_ku(i_burst,1) = bin2dec([bit0 bit15 bit68 bit9 bit10 bit11to15]);
end
l1a_mcd_ku = zeros(N_total_bursts_sar_ku_isp,1); % TBD

%----------D. Altimeter range variables -----------------------------------  
altimeter_range_calibrated_ku = uint32((win_delay_sar_ku_ref * c_cst/2 - 1.3e6) * 1e4);
range_corr_internal_delay_ku = int16(int_delay_cor_cal1 * c_cst/2 * 1e4);
range_corr_com_ku = int16(-mean(z_cog_ant_corr.') * 1e4);

%----------E. Altimeter power variables -----------------------------------
attenuator_calibrated_ku = int16(att_conv_sar_ku * 1e2);
altimeter_power_drift_ku = int16(power_var_cal1_sar_ku * 1e2);
power_corr_digital_processing_ku = int16((onboard_proc_sar)*1e2);
power_scaling_to_antenna_ku = int16(gain_corr_instr_sar * 1e2);

%----------F. Altimeter engineering variables -----------------------------
altimeter_clock_ku = int32(1./T0_sar - 3.95e8)* 1e9;
tm_h0_ku = uint32(h0_sar_isp);
tm_cor2_ku = int16(cor2_sar_isp);
cai_ku = uint32(cai_namb_sar(:,1)).';
fai_ku = uint32(fai_sar(:,1)).';
tm_pri_ku = uint32(pri_sar_isp);
tm_ambiguity_rank_ku = uint16(ambiguity_order_sar_isp);
tm_nimp_ku = uint16(nimp_sar_isp);

%----------L. Waveform related variables -----------------------------------
tm_burst_num_ku = uint8(burst_sar_isp);
i_wfm_cal_corrected = real(wfm_cal_gain_corrected);
q_wfm_cal_corrected = imag(wfm_cal_gain_corrected);
i_scale_factor_ku = zeros(N_ku_pulses_burst_chd,N_total_bursts_sar_ku_isp);
q_scale_factor_ku = zeros(N_ku_pulses_burst_chd,N_total_bursts_sar_ku_isp);
i_samples_ku = int8(zeros(N_samples_sar_chd,N_ku_pulses_burst_chd,N_total_bursts_sar_ku_isp));
q_samples_ku = int8(zeros(N_samples_sar_chd,N_ku_pulses_burst_chd,N_total_bursts_sar_ku_isp));


for i_burst=1:N_total_bursts_sar_ku_isp
    progressbar([],[],[],[],[],[],[],[],[],[],[],[],[],[],(i_burst)/N_total_bursts_sar_ku_isp);
    for i_pulse=1:N_ku_pulses_burst_chd
        i_scale_factor_ku(i_pulse,i_burst) = single(max(abs(i_wfm_cal_corrected(i_burst,i_pulse,:)))*10^0.3 / (2^7-1));
        q_scale_factor_ku(i_pulse,i_burst) = single(max(abs(q_wfm_cal_corrected(i_burst,i_pulse,:)))*10^0.3 / (2^7-1));
        i_samples_ku(:,i_pulse,i_burst) = int8(round(i_wfm_cal_corrected(i_burst,i_pulse,:)*10^0.3./i_scale_factor_ku(i_pulse,i_burst)));
        q_samples_ku(:,i_pulse,i_burst) = int8(round(q_wfm_cal_corrected(i_burst,i_pulse,:)*10^0.3./q_scale_factor_ku(i_pulse,i_burst)));
    end
end

snr_estimation_ku = zeros(N_total_bursts_sar_ku_isp,1); 

%% PACKING L1A

%----------A. Time variables ----------------------------------------------



l1_mode_id_ku_name = 'l1_mode_id_ku';
nccreate(filename,l1_mode_id_ku_name,format_key,netcdf_v4_format,data_type_key,uint8_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite( filename,l1_mode_id_ku_name,l1_mode_id_ku);
ncwriteatt(filename,['/' l1_mode_id_ku_name],long_name_att,'L1 mode ID (Ku-band)');
ncwriteatt(filename,['/' l1_mode_id_ku_name],comment_att,'L1 Mode Identifier. Each L1 record type has a unique ID, as follows: 0 = Null record, 1 = LRM L1b Ku, 2 = LRM L1b C, 3 = SAR Raw L1A Ku, 4 = SAR Raw L1B-S Ku, 5 = SAR Raw L1B Ku, 6 = SAR RMC L1A Ku, 7 = SAR RMC L1B-S Ku, 8 = SAR RMC L1B Ku, 9 = CAL1 LRM L1b Ku, 10 = CAL1 LRM L1b C, 11 = CAL1 SAR L1b Ku, 12 = CAL1 SAR L1b C, 13 = CAL1 INSTR L1b Ku, 14 = CAL1 INSTR L1b C, 15 = Autocal L1b, 16 = CAL2 L1b Ku, 17 = CAL2 L1b C');

time_day_ku_name = 'time_day_ku';
nccreate(filename,time_day_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite( filename,time_day_ku_name,time_day_ku);
ncwriteatt(filename,['/' time_day_ku_name],long_name_att,'Day from 1 Jan 2000 (Ku-band)');
ncwriteatt(filename,['/' time_day_ku_name],units_att,day_units);
ncwriteatt(filename,['/' time_day_ku_name],comment_att,'It contains the day from 1 Jan 2000. Time is in TAI. Time refers to the instant the (theoretical) middle of the burst touches the surface.');

time_seconds_ku_name = 'time_seconds_ku';
nccreate(filename,time_seconds_ku_name,format_key,netcdf_v4_format,data_type_key,uint64_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite( filename,time_seconds_ku_name,time_seconds_ku);
ncwriteatt(filename,['/' time_seconds_ku_name],long_name_att,'Seconds of the day, with microsecond resolution (Ku-band)');
ncwriteatt(filename,['/' time_seconds_ku_name],units_att,seconds_units);
ncwriteatt(filename,['/' time_seconds_ku_name],scale_factor_att,'1e-6');
ncwriteatt(filename,['/' time_seconds_ku_name],comment_att,'Second of the day, with a scaling factor of 1.0e-6, providing a microsecond resolution. Time is in TAI. Time refers to the instant the (theoretical) middle of the burst touches the surface.');

tm_source_sequence_counter_ku_name = 'tm_source_sequence_counter_ku';
nccreate(filename,tm_source_sequence_counter_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,tm_source_sequence_counter_ku_name,tm_source_sequence_counter_ku);
ncwriteatt(filename,['/' tm_source_sequence_counter_ku_name],long_name_att,'Instrument source sequence counter (Ku-band)');
ncwriteatt(filename,['/' tm_source_sequence_counter_ku_name],units_att,number_units);
ncwriteatt(filename,['/' tm_source_sequence_counter_ku_name],comment_att,'This is the instrument record counter, copied from ISP. Each L1A SAR record (dimension ku_rec) corresponds to 1 single ISP (SAR Raw or SAR RMC');

burst_counter_ku_name = 'burst_counter_ku';
nccreate(filename,burst_counter_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,burst_counter_ku_name,burst_counter_ku);
ncwriteatt(filename,['/' burst_counter_ku_name],long_name_att,'L1A record burst counter (Ku-band)');
ncwriteatt(filename,['/' burst_counter_ku_name],units_att,number_units);
ncwriteatt(filename,['/' burst_counter_ku_name],comment_att,' The L1A record counter, starting from 0 for each product');

%----------B. Orbit and attitude variables ------------------------------
latitude_ku_name = 'latitude_ku';
nccreate(filename,latitude_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,latitude_ku_name,latitude_ku);
ncwriteatt(filename,['/' latitude_ku_name],long_name_att,'latitude (positive N, negative S) (Ku-band)');
ncwriteatt(filename,['/' latitude_ku_name],units_att,degrees_units);
ncwriteatt(filename,['/' latitude_ku_name],scale_factor_att,'1e-7');
ncwriteatt(filename,['/' latitude_ku_name],comment_att,'Latitude of measurement [-90, +90]: Positive at Nord, Negative at South');

longitude_ku_name = 'longitude_ku';
nccreate(filename,longitude_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,longitude_ku_name,longitude_ku);
ncwriteatt(filename,['/' longitude_ku_name],long_name_att,'Longitude (degrees East) (Ku-band)');
ncwriteatt(filename,['/' longitude_ku_name],units_att,degrees_units);
ncwriteatt(filename,['/' longitude_ku_name],scale_factor_att,'1e-7');
ncwriteatt(filename,['/' longitude_ku_name],comment_att,'Longitude of measurement [0, 360]: Positive at East');

com_altitude_ku_name = 'com_altitude_ku';
nccreate(filename,com_altitude_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,com_altitude_ku_name,com_altitude_ku);
ncwriteatt(filename,['/' com_altitude_ku_name],long_name_att,'CoM altitude (Ku-band)');
ncwriteatt(filename,['/' com_altitude_ku_name],units_att,meters_units);
ncwriteatt(filename,['/' com_altitude_ku_name],add_offset_att,'1.3e6');
ncwriteatt(filename,['/' com_altitude_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' com_altitude_ku_name],comment_att,'Altitude of the satellite Centre of Mass.');

com_altitude_rate_ku_name = 'com_altitude_rate_ku';
nccreate(filename,com_altitude_rate_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,com_altitude_rate_ku_name,com_altitude_rate_ku);
ncwriteatt(filename,['/' com_altitude_rate_ku_name],long_name_att,'CoM altitude rate (Ku-band)');
ncwriteatt(filename,['/' com_altitude_rate_ku_name],units_att,meters_per_second_units);
ncwriteatt(filename,['/' com_altitude_rate_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' com_altitude_rate_ku_name],comment_att,'Instantaneous altitude rate at the Centre of Mass');

com_velocity_vector_ku_name = 'com_velocity_vector_ku';
nccreate(filename,com_velocity_vector_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{space_3D_dimension ,space_3D_dimension_size,ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,com_velocity_vector_ku_name,com_velocity_vector_ku);
ncwriteatt(filename,['/' com_velocity_vector_ku_name],long_name_att,'Velocity vector in ITRF, components: [1] x, [2] y, [3] z (Ku-band)');
ncwriteatt(filename,['/' com_velocity_vector_ku_name],units_att,meters_per_second_units);
ncwriteatt(filename,['/' com_velocity_vector_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' com_velocity_vector_ku_name],comment_att,'Velocity Vector at the centre of mass in ITRF. The 3 components are given according to the "space_3D" dimension: [1] x, [2] y, [3] z');

satellite_mispointing_ku_name = 'satellite_mispointing_ku';
nccreate(filename,satellite_mispointing_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{space_3D_dimension ,space_3D_dimension_size,ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,satellite_mispointing_ku_name,satellite_mispointing_ku);
ncwriteatt(filename,['/' satellite_mispointing_ku_name],long_name_att,'Mispointing angle, measures by STRs: [1] Roll, [2] Pitch, [3] Yaw (Ku-band)');
ncwriteatt(filename,['/' satellite_mispointing_ku_name],units_att,degrees_units);
ncwriteatt(filename,['/' satellite_mispointing_ku_name],scale_factor_att,'1e-7');
ncwriteatt(filename,['/' satellite_mispointing_ku_name],comment_att,'Attitude mispointing, measured by STRs and post-processed by AOCS or by ground facility. The 3 components are given according to the "space_3D" dimension: [1] Roll, [2] Pitch, [3] Yaw. This variable includes the "mispointing bias" given by the variable mispointing_bias_ku. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: roll = pitch = yaw = 0');

mispointing_bias_ku_name = 'mispointing_bias_ku';
nccreate(filename,mispointing_bias_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{space_3D_dimension, space_3D_dimension_size});
ncwrite(filename,mispointing_bias_ku_name,mispointing_bias_ku);
ncwriteatt(filename,['/' mispointing_bias_ku_name],long_name_att,'Mispointing bias from in-flight calibration: [1] Roll, [2] Pitch, [3] Yaw (Ku-band)');
ncwriteatt(filename,['/' mispointing_bias_ku_name],units_att,degrees_units);
ncwriteatt(filename,['/' mispointing_bias_ku_name],scale_factor_att,'1e-7');
ncwriteatt(filename,['/' mispointing_bias_ku_name],comment_att,'The mispointing bias will correct for structural misalignment not accounted for during on-ground characterization. It will provide absolute mispointing value. It will be estimated using altimeter''s measurements data. The mispointing bias is given for [1] roll, [2] pitch and [3] yaw and it is constant within the product.');

%----------C. Configuration and quality variables -------------------------
l1_instrument_configuration_ku_name = 'l1_instrument_configuration_ku';
nccreate(filename,l1_instrument_configuration_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,l1_instrument_configuration_ku_name,l1_instrument_configuration_ku);
ncwriteatt(filename,['/' l1_instrument_configuration_ku_name],long_name_att,'Instrument configuration L1 (Ku-band)');
ncwriteatt(filename,['/' l1_instrument_configuration_ku_name],comment_att,'It contains flags from ISPs related to P4 configuration and tracking. For more details, see ESA IODD "JC-ID-ESA-GP-0059"');

l1a_mcd_ku_name = 'l1a_mcd_ku';
nccreate(filename,l1a_mcd_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
% ncwrite(filename,l1a_mcd_ku_name,l1a_mcd_ku);
ncwriteatt(filename,['/' l1a_mcd_ku_name],long_name_att,'Measurement confidence data (MCD) L1A (Kuband)');
ncwriteatt(filename,['/' l1a_mcd_ku_name],comment_att,'It contains flags related to the L1A processing chain. For more details, see ESA IODD "JC-ID-ESA-GP-0059"');

%----------D. Altimeter range variables -----------------------------------
altimeter_range_calibrated_ku_name = 'altimeter_range_calibrated_ku';
nccreate(filename,altimeter_range_calibrated_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,altimeter_range_calibrated_ku_name,altimeter_range_calibrated_ku);
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],long_name_att,'Calibrated 1-way range: CoM to middle range window (at sample Ns/2 from 0) (Ku-band)');
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],units_att,meters_units);
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],add_offset_att,'1.3e6');
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' altimeter_range_calibrated_ku_name],comment_att,'This is the 1-way distance from the satellite''s Center of Mass to the middle of the range window (sample Ns/2 from 0). It includes the following range calibrations: (a) range_corr_internal_delay_ku, (b) range_corr_com_ku. Note: the actual altimeter clock (variable "altimeter_clock_ku") has been used to compute the altimeter range');

range_corr_internal_delay_ku_name = 'range_corr_internal_delay_ku';
nccreate(filename,range_corr_internal_delay_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,range_corr_internal_delay_ku_name,range_corr_internal_delay_ku);
ncwriteatt(filename,['/' range_corr_internal_delay_ku_name],long_name_att,'Range correction:internal path delay (waveguides + CAL1) (Ku-band)');
ncwriteatt(filename,['/' range_corr_internal_delay_ku_name],units_att,meters_units);
ncwriteatt(filename,['/' range_corr_internal_delay_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' range_corr_internal_delay_ku_name],comment_att,'The internal path delay is the sum of 2 components: (a) waveguides component: constant part measured pre-launch, (b) CAL1 P4 internal delay: periodically measured onflight');

range_corr_com_ku_name = 'range_corr_com_ku';
nccreate(filename,range_corr_com_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,range_corr_com_ku_name,range_corr_com_ku);
ncwriteatt(filename,['/' range_corr_com_ku_name],long_name_att,'Range correction: distance CoM antenna(Ku-band)');
ncwriteatt(filename,['/' range_corr_com_ku_name],units_att,meters_units);
ncwriteatt(filename,['/' range_corr_com_ku_name],scale_factor_att,'1e-4');
ncwriteatt(filename,['/' range_corr_com_ku_name],comment_att,'It is the distance in the z-component between the Centre of Mass and the Antenna reference point. The Centre of Mass is assumed constant in the GPP, although it is likely to vary in P4');

%----------E. Altimeter power variables -----------------------------------
attenuator_calibrated_ku_name = 'attenuator_calibrated_ku';
nccreate(filename,attenuator_calibrated_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,attenuator_calibrated_ku_name,attenuator_calibrated_ku);
ncwriteatt(filename,['/' attenuator_calibrated_ku_name],long_name_att,'Power scaling: ATT calibrated (from ATT table) (Ku-band)');
ncwriteatt(filename,['/' attenuator_calibrated_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' attenuator_calibrated_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' attenuator_calibrated_ku_name],comment_att,'The ISPs contain the ATT command applied within a radar cycle. The actual ATT value, that slightly differs from the nominal one, is extracted from a look-up table');

altimeter_power_drift_ku_name = 'altimeter_power_drift_ku';
nccreate(filename,altimeter_power_drift_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,altimeter_power_drift_ku_name,altimeter_power_drift_ku);
ncwriteatt(filename,['/' altimeter_power_drift_ku_name],long_name_att,'Power scaling: power instrument aging correction (CAL1) (Ku-band)');
ncwriteatt(filename,['/' altimeter_power_drift_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' altimeter_power_drift_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' altimeter_power_drift_ku_name],comment_att,'This is a measure of the instrument power decay over the time, regularly measured from CAL1 impulse response with respect to a reference at the beginning of the mission.');

power_corr_digital_processing_ku_name = 'power_corr_digital_processing_ku';
nccreate(filename,power_corr_digital_processing_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type);
ncwrite(filename,power_corr_digital_processing_ku_name,power_corr_digital_processing_ku);
ncwriteatt(filename,['/' power_corr_digital_processing_ku_name],long_name_att,'Power scaling: onboard digital processing (Ku-band)');
ncwriteatt(filename,['/' power_corr_digital_processing_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' power_corr_digital_processing_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' power_corr_digital_processing_ku_name],comment_att,'Each measurement mode (LRM, SAR, RMC) is subject to a specific on-board digital processing that affects the signal scaling. This factor includes all the digital processing scaling factors, which are constant for each mode.');

power_scaling_to_antenna_ku_name = 'power_scaling_to_antenna_ku';
nccreate(filename,power_scaling_to_antenna_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,power_scaling_to_antenna_ku_name,power_scaling_to_antenna_ku);
ncwriteatt(filename,['/' power_scaling_to_antenna_ku_name],long_name_att,'Power scaling: overall gain from raw waveform to antenna flange (Ku-band)');
ncwriteatt(filename,['/' power_scaling_to_antenna_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' power_scaling_to_antenna_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' power_scaling_to_antenna_ku_name],comment_att,'The overall instrument gain, including: (a) attenuator_calibrated_ku, (b) altimeter_power_drift_ku, (c) power_corr_digital_processing_ku, (d) RFU gains. Antenna gain is excluded. This variable is applied to the L1A waveforms as power scaling: the resulting power is the received power at the antenna flange.');

%----------F. Altimeter engineering variables -----------------------------
altimeter_clock_ku_name = 'altimeter_clock_ku';
nccreate(filename,altimeter_clock_ku_name,format_key,netcdf_v4_format,data_type_key,int32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,altimeter_clock_ku_name,altimeter_clock_ku);
ncwriteatt(filename,['/' altimeter_clock_ku_name],long_name_att,'Altimeter clock (Ku-band)');
ncwriteatt(filename,['/' altimeter_clock_ku_name],units_att,Hz_units);
ncwriteatt(filename,['/' altimeter_clock_ku_name],add_offset_att,'395e6');
ncwriteatt(filename,['/' altimeter_clock_ku_name],scale_factor_att,'1e-9');
ncwriteatt(filename,['/' altimeter_clock_ku_name],comment_att,'This is the actual altimeter clock. The altimeter clock is based upon the USO clock. The nominal USO clock is 10MHz, while the nominal altimeter clock is 395 MHz. The actual USO clock is provided regularly as a drift information: based on the USO drift, the actual altimeter clock is computed, assuming a linear dependency.');

tm_h0_ku_name = 'tm_h0_ku';
nccreate(filename,tm_h0_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,tm_h0_ku_name,tm_h0_ku);
ncwriteatt(filename,['/' tm_h0_ku_name],long_name_att,'H0 telemetry (same in 1 radar cycle) (Ku-band)');
ncwriteatt(filename,['/' tm_h0_ku_name],units_att,T0d64_units);
ncwriteatt(filename,['/' tm_h0_ku_name],comment_att,'The H0 (initial altitude instruction) value, copied from ISP. The altimeter provides one H0 value per radar cycle. It is used in combination with the COR2 to compute the cai_ku and fai_ku instructions. NOTE: T0 = 1/altimeter_clock_ku_ku.');

tm_cor2_ku_name = 'tm_cor2_ku';
nccreate(filename,tm_cor2_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,tm_cor2_ku_name,tm_cor2_ku);
ncwriteatt(filename,['/' tm_cor2_ku_name],long_name_att,'COR2 telemetry (same in 1 radar cycle) (Ku-band)');
ncwriteatt(filename,['/' tm_cor2_ku_name],units_att,T0d16d64_units);
ncwriteatt(filename,['/' tm_cor2_ku_name],comment_att,'The COR2 (altitude rate estimation) value, copied from ISP. The altimeter estimates one COR2 value per radar cycle. It is used in combination with the H0 to compute the cai_ku and fai_ku instructions. NOTE: T0 = 1/altimeter_clock_ku_ku.');

cai_ku_name = 'cai_ku';
nccreate(filename,cai_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,cai_ku_name,cai_ku);
ncwriteatt(filename,['/' cai_ku_name],long_name_att,'Coarse Altitude Instruction (Same for 64 pulses)');
ncwriteatt(filename,['/' cai_ku_name],units_att,T0x4_units);
ncwriteatt(filename,['/' cai_ku_name],comment_att,'This is the Coarse Altitude Instruction of the burst. The cai_ku can potentially change within each pulse within the burst. The L1A processor has the (configurable) capability of removing the cai_ku variation within the burst and replace it with the altitude rate estimated from the orbit. The L1A product gives the cai_ku of the first pulse within the burst. NOTE: T0 = 1/altimeter_clock_ku_ku');

fai_ku_name = 'fai_ku';
nccreate(filename,fai_ku_name,format_key,netcdf_v4_format,data_type_key,uint32_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,fai_ku_name,fai_ku);
ncwriteatt(filename,['/' fai_ku_name],long_name_att,'Fine Altitude Instruction (Same for 64 pulses)');
ncwriteatt(filename,['/' fai_ku_name],units_att,T0d16d64_units);
ncwriteatt(filename,['/' fai_ku_name],comment_att,'This is the Fine Altitude Instruction of the burst. The fai_ku can potentially change within each pulse within the burst. The L1A processor has the (configurable) capability of removing the fai_ku variation within the burst and replace it with the altitude rate estimated from the orbit. The L1A product gives the fai_ku of the first pulse within the burst. NOTE: T0 = 1/altimeter_clock_ku_ku');

tm_pri_ku_name = 'tm_pri_ku';
nccreate(filename,tm_pri_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,tm_pri_ku_name,tm_pri_ku);
ncwriteatt(filename,['/' tm_pri_ku_name],long_name_att,'Pulse Repetition Interval (form ISP)');
ncwriteatt(filename,['/' tm_pri_ku_name],units_att,T0x8_units);
ncwriteatt(filename,['/' tm_pri_ku_name],comment_att,'The "Pulse Repetition Interval. PRI is constant within all received pulses in a radar cycle, but it can change within consecutive radar cycles. It is provided in counters of (T0*8) by the altimeter, with T0 = 1/altimeter_clock_ku_ku.');

tm_ambiguity_rank_ku_name = 'tm_ambiguity_rank_ku';
nccreate(filename,tm_ambiguity_rank_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,tm_ambiguity_rank_ku_name,tm_ambiguity_rank_ku);
ncwriteatt(filename,['/' tm_ambiguity_rank_ku_name],long_name_att,'Ambiguity Rank (from ISP)');
ncwriteatt(filename,['/' tm_ambiguity_rank_ku_name],units_att,number_units);
ncwriteatt(filename,['/' tm_ambiguity_rank_ku_name],comment_att,'The Pulse Ambiguity Rank, taken from ISP: it is the number of pulses that are transmitted between the transmission and the reception of each pulse.');

tm_nimp_ku_name = 'tm_nimp_ku';
nccreate(filename,tm_nimp_ku_name,format_key,netcdf_v4_format,data_type_key,uint16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,tm_nimp_ku_name,tm_nimp_ku);
ncwriteatt(filename,['/' tm_nimp_ku_name],long_name_att,'Number of pulses in 1 radar cycle (from ISP)');
ncwriteatt(filename,['/' tm_nimp_ku_name],units_att,number_units);
ncwriteatt(filename,['/' tm_nimp_ku_name],comment_att,'It is the number of pulses in 1 radar cycle, taken from ISP');

%----------L. Waveform related variables -----------------------------
tm_burst_num_ku_name = 'tm_burst_num_ku';
nccreate(filename,tm_burst_num_ku_name,format_key,netcdf_v4_format,data_type_key,uint8_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,tm_burst_num_ku_name,tm_burst_num_ku);
ncwriteatt(filename,['/' tm_burst_num_ku_name],long_name_att,'burst_number (from 1 to 7)');
ncwriteatt(filename,['/' tm_burst_num_ku_name],units_att,number_units);
ncwriteatt(filename,['/' tm_burst_num_ku_name],comment_att,'Each radar cycle contains 7 bursts: this is the burst counter within the radar cycle, from 1 to 7. NOTE: the L1A contains full radar cycles(groups of 7 bursts), which makes the dimension ku_rec multiple of 7.');

i_samples_ku_name = 'i_samples_ku';
nccreate(filename,i_samples_ku_name,format_key,netcdf_v4_format,data_type_key,int8_type,dimensions_key,{ns_dimension,N_samples_sar_chd,np_dimension,N_ku_pulses_burst_chd,ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,i_samples_ku_name,i_samples_ku);
ncwriteatt(filename,['/' i_samples_ku_name],long_name_att,'SAR L1A Ku i-samples, arranged in bursts of 64x256 elements. I samples are scaled to range [-127,+127].');
ncwriteatt(filename,['/' i_samples_ku_name],units_att,number_units);
ncwriteatt(filename,['/' i_samples_ku_name],comment_att,'The i component of each L1A waveform. Each waveform is a fully calibrated, pulse limited complex echo. Echoes are phase coherent, which makes the SAR processing possible. Each pulse within the burst is: (a) given in the time domain, (b) power calibrated with "power_scaling_to_antenna_ku", (c) calibrated with CAL2 mask (when correction switched on). All pulses within a burst are: (d) aligned in range, with the altitude rate being removed by the altimeter or by the L1A Processor, (e) calibrated with CAL1 pulse-to-pulse amplitude and phase correction (when correction switched on). A final scaling, given in the variable "i_scale_factor_ku", is applied in order to best fit the i component into 1 byte.');

q_samples_ku_name = 'q_samples_ku';
nccreate(filename,q_samples_ku_name,format_key,netcdf_v4_format,data_type_key,int8_type,dimensions_key,{ns_dimension,N_samples_sar_chd,np_dimension,N_ku_pulses_burst_chd,ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,q_samples_ku_name,q_samples_ku);
ncwriteatt(filename,['/' q_samples_ku_name],long_name_att,'SAR L1A Ku q-samples, arranged in bursts of 64x256 elements. Q samples are scaled to range [-127,+127].');
ncwriteatt(filename,['/' q_samples_ku_name],units_att,number_units);
ncwriteatt(filename,['/' q_samples_ku_name],comment_att,'The q component of each L1A waveform. Each waveform is a fully calibrated, pulse limited complex echo. Echoes are phase coherent, which makes the SAR processing possible. Each pulse within the burst is: (a) given in the time domain, (b) power calibrated with "power_scaling_to_antenna_ku", (c) calibrated with CAL2 mask (when correction switched on). All pulses within a burst are: (d) aligned in range, with the altitude rate being removed by the altimeter or by the L1A Processor, (e) calibrated with CAL1 pulse-to-pulse amplitude and phase correction (when correction switched on). A final scaling, given in the variable "q_scale_factor_ku", is applied in order to best fit the q component into 1 byte.');

i_scale_factor_ku_name = 'i_scale_factor_ku';
nccreate(filename,i_scale_factor_ku_name,format_key,netcdf_v4_format,data_type_key,float_type,dimensions_key,{np_dimension,N_ku_pulses_burst_chd,ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,i_scale_factor_ku_name,i_scale_factor_ku);
ncwriteatt(filename,['/' i_scale_factor_ku_name],long_name_att,'I scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
ncwriteatt(filename,['/' i_scale_factor_ku_name],units_att,sqrtW_per_count_units);
ncwriteatt(filename,['/' i_scale_factor_ku_name],comment_att,'The i_samples_ku scaling factor, computed in order to best fit the i_samples_ku within 1 byte. The scaling, needed to convert the i_samples_ku into sqrt(watt), is applied as follows: i_samples_ku_sqr_watt(ku_rec,Np, Ns) = i_samples_ku(ku_rec,Np, Ns) * i_scale_factor_ku(ku_rec,Np) ');

q_scale_factor_ku_name = 'q_scale_factor_ku';
nccreate(filename,q_scale_factor_ku_name,format_key,netcdf_v4_format,data_type_key,float_type,dimensions_key,{np_dimension,N_ku_pulses_burst_chd,ku_rec_dimension, N_total_bursts_sar_ku_isp});
ncwrite(filename,q_scale_factor_ku_name,q_scale_factor_ku);
ncwriteatt(filename,['/' q_scale_factor_ku_name],long_name_att,'Q scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
ncwriteatt(filename,['/' q_scale_factor_ku_name],units_att,sqrtW_per_count_units);
ncwriteatt(filename,['/' q_scale_factor_ku_name],comment_att,'The q_samples_ku scaling factor, computed in order to best fit the q_samples_ku within 1 byte. The scaling, needed to convert the q_samples_ku into sqrt(watt), is applied as follows: q_samples_ku_sqr_watt(ku_rec,Np, Ns) = q_samples_ku(ku_rec,Np, Ns) * q_scale_factor_ku(ku_rec,Np)');

snr_estimation_ku_name = 'snr_estimation_ku';
nccreate(filename,snr_estimation_ku_name,format_key,netcdf_v4_format,data_type_key,int16_type,dimensions_key,{ku_rec_dimension, N_total_bursts_sar_ku_isp});
% ncwrite(filename,snr_estimation_ku_name,snr_estimation_ku);
ncwriteatt(filename,['/' snr_estimation_ku_name],long_name_att,'SNR estimation on L1A burst');
ncwriteatt(filename,['/' snr_estimation_ku_name],units_att,dB_units);
ncwriteatt(filename,['/' snr_estimation_ku_name],scale_factor_att,'0.01');
ncwriteatt(filename,['/' snr_estimation_ku_name],comment_att,'SNR estimation done on L1A burst-averaged waveform as: peak / noise_level, where noise_level is estimated from the first n (TBC) samples');

ncwriteatt(filename,'/','creation_time',datestr(now));
ncwriteatt(filename,'/','data_info','data simulated by ESA and processed by isardSAT');


time = toc(t6);
minutes_reading = floor(time/60);
secs_reading = time - minutes_reading*60;
disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1A']);

