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
% READFBRrecord: function that reads the FBR data set from input filename,
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
%   L1A record           
% ----------------------------------------------------------
% 
% Author:   Albert Garcia / isardSAT
%
% Version  record
% 1.0 2016/03/22 Cloned from readFBR_CryoSat_SAR
% 2.0 2016/04/28 Open the binary file just once and use the position of the
% last pointer
% 2.1 2016/05/18 Removed ds_offset and dsr_size as globals
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [burst] = readFBR_CryoSat_SAR_record(files,burst, i_burst)

global mode N_samples 
global flat_coeff_cst semi_major_axis_cst
global pri_sar_chd bri_chd prf_sar_chd brf_chd T0_chd
global pri_sar_nom bri_sar_nom prf_sar_nom brf_sar_nom c_cst
global height_rate_application_cnf FAI_application_cnf

t1 = tic;
switch mode
    case 'SAR'
        % Fields mapping size            51               52  53
        waveform_group_record_size     = 64*N_samples*2+  2+  2;
    case 'SIN'
       % Fields mapping size            54               55               56  57
        waveform_group_record_size     = 64*N_samples*2+  64*N_samples*2+  2+  2;
end
%find the memory positions of the record to be read

%% From Product Specs
% Fields mapping size           1   2 	3   4   5   6   7   8   9   10  11      12      13      14
time_orbit_group_record_size =  12+ 4+  2+  2+  4+  4+  4+  4+  4+  4+  3*4+    3*4+    3*4+    4;
% Fields mapping size           15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31	32 	33	34
measurements_group_record_size = 8+ 4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4;
% Fields mapping size           35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50
% corrections_group_record_size  = 4+ 4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4+  4;

i_burst_db = ceil(i_burst/20);
i_burst_record = mod(i_burst-1,20);
remaining_records = (20-(i_burst_record+1)); 
burst.source_seq_count_sar_isp=i_burst_db;


%% READ FBR record        
fseek(files.fid,files.sph.dsds(1).ds_offset+(i_burst_db-1)*files.sph.dsds(1).dsr_size,'bof'); % last records of the product


%% READ TIME ORBIT
fseek(files.fid,(i_burst_record)*time_orbit_group_record_size,'cof');
[burst] = read_FBR_CryoSat_time_orbit_group(files.fid,burst);
fseek(files.fid,(remaining_records)*time_orbit_group_record_size,'cof');

%% READ MEASUREMENTS 
fseek(files.fid,(i_burst_record)*measurements_group_record_size,'cof');
[burst] = read_FBR_CryoSat_measurements_group(files.fid, burst);
fseek(files.fid,(remaining_records)*measurements_group_record_size,'cof');

burst.T0_sar     = T0_chd*(burst.USO_correction+1); %T0 real
pri_sar_nom = pri_sar_chd;
burst.pri_sar = pri_sar_nom;
prf_sar_nom = prf_sar_chd; 
bri_sar_nom = bri_chd; 
brf_sar_nom = brf_chd; 





%% READ GEOPHYSICAL CORRECTTIONS
[burst] = read_FBR_CryoSat_corrections_group(files.fid, burst);

%% READ WAVEFORMS
fseek(files.fid,(i_burst_record)*waveform_group_record_size,'cof');
[burst] = read_FBR_CryoSat_waveform_group(files.fid, burst);

    
%% ALIGN with S6 L1A
burst.USO_correction    = burst.win_delay_sar_ku.*(burst.USO_correction); % seconds correction for window delay
burst.win_delay_sar_ku  = burst.win_delay_sar_ku + burst.USO_correction + burst.instrument_range_correction_tx_rx /(c_cst/2); % seconds
if strcmp(mode,'SIN')
    burst.win_delay_sar_ku_2 = burst.win_delay_sar_ku-burst.instrument_range_correction_tx_rx /(c_cst/2) + burst.instrument_range_correction_rx /(c_cst/2); % seconds
end

%% Perform co-registration of two channels for SARin
if strcmp(mode,'SIN')
    [burst]=coregistration_channels_SARin(burst);
end



burst.time_sar_ku = burst.days + burst.seconds + burst.microseconds;

p = lla2ecef([burst.lat_sar_sat.',burst.lon_sar_sat.',burst.alt_sar_sat.'],flat_coeff_cst,semi_major_axis_cst);
burst.x_sar_sat = p(:,1).';
burst.y_sar_sat = p(:,2).';
burst.z_sar_sat = p(:,3).';

burst.lat_sar_surf = burst.lat_sar_sat;
burst.lon_sar_surf = burst.lon_sar_sat;
burst.alt_sar_surf = burst.alt_sar_sat - burst.win_delay_sar_ku * c_cst/2;


% geod2cart(SURF)
p = lla2ecef([burst.lat_sar_surf,burst.lon_sar_surf,burst.alt_sar_surf],flat_coeff_cst,semi_major_axis_cst);
burst.x_sar_surf = p(1).';
burst.y_sar_surf = p(2).';
burst.z_sar_surf = p(3).';

[~,burst.doppler_ang_sar_sat] = compute_height_rate(1, burst.x_vel_sat_sar, burst.y_vel_sat_sar, burst.z_vel_sat_sar,burst.x_sar_sat ,burst.y_sar_sat ,burst.z_sar_sat ,burst.x_sar_surf,burst.y_sar_surf,burst.z_sar_surf);






%temps = toc(t1);
%minuts = floor(temps/60);
%segs = temps - minuts*60;
% disp(['Han passat ',num2str(minuts),' minuts i ',num2str(segs),' segons']);

end