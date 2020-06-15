%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-1BS product using Sentinel-3 like format
% Ref: S3-TN-ESA-SR-0433 SRAL L1A-1BS IODD V1.4
%
% ---------------------------------------------------------
% Objective: Pack variables and write the NETCDF
% 
% INPUTs : Workspace
% OUTPUTs: TM Structure as defined on isardSAT_JasonCS_DPM
%
% ----------------------------------------------------------
% Author:    Eduard Makhoul/ isardSAT
%            Gorka Moyano  / isardSAT
%            Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT ()
% 
% Versions
% 1.0 
% 1.1 Updated time conversion for data between 2010 and 2016 (CR2)
% 2.0 Transformed to a function. Writting one record
% 2.1 case 'SIN' added, changed N_samples_sar_chd by N_samples
% 2.2 added cnf flag processing_mode_cnf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function [files] = create_NetCDF_L1BS_S3(files, N_surfs_loc_estimated)

global N_samples   semi_major_axis_cst flat_coeff_cst;
global bw_ku_chd zp_fact_range_cnf  sec_in_day_cst pi_cst
global mission mode N_max_beams_stack_chd processing_mode_cnf

% t6 = tic;

date_creation = datestr(now, '_yyyymmddTHHMMSS_');
switch mission
    case 'CR2'        
        switch mode
            case 'SAR'                
                files.filename_netCDF_BS = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA_BS_',...
                                files.sph.product_info.product_id(20:20+30),...
                                date_creation,...
                                'isd','.nc');
            case 'SIN'                
                files.filename_netCDF_BS = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA_BS_',...
                                files.sph.product_info.product_id(20:20+30),...
                                date_creation,...
                                'isd','.nc');
        end
    case {'S3','S3_','S3A','S3B'}
        files.filename_netCDF_BS = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA_BS_',...
                                date_creation,...
                                'isd','.nc');
    case {'S6_'} 
end
ncid = netcdf.create(files.filename_netCDF_BS,'CLOBBER');

long_name_att = 'long_name';
std_name_att = 'standard_name';
calendar_name_att='calendar';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';
flag_values_att='flag_values';
flag_desc_att='flag_meanings';
dimensions_key = 'Dimensions';
format_key = 'Format';
data_type_key = 'DataType';
fill_value_key='FillValue';

netcdf_v4_format = 'netcdf4';

% ku_rec_dimension = 'time_l1b_echo_sar_ku';
% nl_dimension = 'Nl';
% nl_dimension_size = N_max_beams_stack;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;

ku_rec_dimension = netcdf.defDim(ncid,'time_l1b_echo_sar_ku',N_surfs_loc_estimated);
nl_dimension = netcdf.defDim(ncid,'Nl',N_max_beams_stack_chd);
ns_dimension = netcdf.defDim(ncid,'Ns',N_samples*zp_fact_range_cnf);
space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);

day_units = 'day';
seconds_units = 'seconds';
seconds_3dot125d64d1e9_units='3.125/64*1e-9 seconds';
seconds_3dot125d1024d1e9_units='3.125/1024*1e-9 seconds';
number_units = 'count';
degrees_units = 'degrees';
meters_units = 'meters';
meters_per_second_units = 'm/s';
rate_per_second_units='1/s';
dB_units = 'dB';
fft_pow_units='FFT power unit';
Hz_units = 'Hz';
T0d64_units = 'T0/64';
T0d16d64_units = 'T0/16/64';
W_per_count_units = 'Watt/#';
sqrtW_per_count_units = 'sqrt(Watt)/#';
rad_units = 'rad';
percent_units='percent';


int8_type = 'NC_BYTE';
uint8_type = 'NC_BYTE';
int16_type = 'NC_SHORT';
uint16_type = 'NC_SHORT';
int32_type = 'NC_INT';
uint32_type = 'NC_INT';
uint64_type = 'uint64';
float_type = 'NC_FLOAT';
double_type= 'NC_DOUBLE';


%% CODING L1B

%--------- Global Attribtues of the netCDF ---------------------------
switch mission
    case 'CR2'
        %altimeter sensor name
        switch files.sph.ins_info.ins_id
            case 'A'
               altimeter_sensor_name='SIRAL Nominal'; 
            case 'B'
               altimeter_sensor_name='SIRAL Redundant'; 
        end
        altimeter_sensor_name=deblank(files.sph.ins_info.ins_conf);
        
        %gnss sensor name
        if isempty(files.sph.gnss_info) 
            gnss_sensor_name='Not available';
        else 
            gnss_sensor_name=files.sph.gnss_info;
        end
        %doris sensor name
        if isempty(files.sph.gnss_info) 
            doris_sensor_name='Not available';
        else 
            doris_sensor_name=files.sph.gnss_info;
        end
        %Acquisition station
        acq_station_name=files.sph.acq_station;
        % UTC date of the first measurement
        first_meas_time=files.sph.product_info.product_time_info.produc_start_time;
        % UTC date of the last measurement
        last_meas_time=files.sph.product_info.product_time_info.produc_stop_time;
        % Name of the altimeter level 0 data file
        xref_altimeter_level0=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'SIRAL_LEVEL_0',length('SIRAL_LEVEL_0'))).filename;
        % Name of the file containing the DORIS-derived USO frequency
        xref_altimeter_orbit=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'ORBIT',length('ORBIT'))).filename;
        % Name of the file containing the DORIS-derived USO frequency
        xref_doris_USO=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'DORIS',length('DORIS'))).filename;
        % Name of the LTM  file containing the SAR mode CAL1 parameters
        xref_altimeter_ltm_sar_cal1=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'CALIBRATION_TYPE_1',length('CALIBRATION_TYPE_1'))).filename;
        % Name of the LTM  file containing the Ku-band CAL2 parameters
        idx=find(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'CALIBRATION_TYPE_2',length('CALIBRATION_TYPE_2'))~=0);
        if isempty(idx)
           xref_altimeter_ltm_ku_cal2='Not available (not applied in the product)';
        else
           xref_altimeter_ltm_ku_cal2=files.sph.dsds(idx).filename;
        end
        % Name of the LTM file containing the C-band CAL2 parameters         
        xref_altimeter_ltm_c_cal2='Not available for CR2';
        % Name of the altimeter characterisation data file
        xref_altimeter_characterisation=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'IPF_RA_DATABASE',length('IPF_RA_DATABASE'))).filename;
        % Semi-major axis of the reference ellipsoid
        semi_major_ellipsoid_axis=num2str(semi_major_axis_cst,15);
        % Flattening coeffcient of the reference ellipsoid
        ellipsoid_flattening=num2str(flat_coeff_cst,15);    
        % Attributes related to the OrbitalN_surfs_loc_estimatedormation required by Porto
        orbit_phase_code=files.sph.orbit_info.phase_code;
        orbit_cycle_num=((files.sph.orbit_info.cycle_num));
        orbit_REL_Orbit=((files.sph.orbit_info.rel_orbit));
        orbit_ABS_Orbit_Start=((files.sph.orbit_info.ABS_Orbit_Start));
        orbit_Rel_Time_ASC_Node_Start=((files.sph.orbit_info.Rel_Time_ASC_Node_Start));
        orbit_ABS_Orbit_Stop=((files.sph.orbit_info.ABS_Orbit_Stop));
        orbit_Rel_Time_ASC_Node_Stop=((files.sph.orbit_info.Rel_Time_ASC_Node_Stop));
        orbit_Equator_Cross_Time=(files.sph.orbit_info.Equator_Cross_Time);
        orbit_Equator_Cross_Long=((files.sph.orbit_info.Equator_Cross_Long));
        orbit_Ascending_Flag=files.sph.orbit_info.Ascending_Flag;
    case {'S3','S3_','S3A','S3B'}
                altimeter_sensor_name='SRAL '; 
        %gnss sensor name
%         if isempty(files.sph.gnss_info) 
%             gnss_sensor_name='Not available';
%         else 
%             gnss_sensor_name=files.sph.gnss_info;
%         end
%         %doris sensor name
%         if isempty(files.sph.gnss_info) 
%             doris_sensor_name='Not available';
%         else 
%             doris_sensor_name=files.sph.gnss_info;
%         end
        %Acquisition station
%         acq_station_name=files.sph.acq_station;
%         % UTC date of the first measurement
%         first_meas_time=files.sph.product_info.product_time_info.produc_start_time;
%         % UTC date of the last measurement
%         last_meas_time=files.sph.product_info.product_time_info.produc_stop_time;
%         % Name of the altimeter level 0 data file
%         xref_altimeter_level0=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'SIRAL_LEVEL_0',length('SIRAL_LEVEL_0'))).filename;
%         % Name of the file containing the DORIS-derived USO frequency
%         xref_altimeter_orbit=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'ORBIT',length('ORBIT'))).filename;
%         % Name of the file containing the DORIS-derived USO frequency
%         xref_doris_USO=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'DORIS',length('DORIS'))).filename;
%         % Name of the LTM  file containing the SAR mode CAL1 parameters
%         xref_altimeter_ltm_sar_cal1=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'CALIBRATION_TYPE_1',length('CALIBRATION_TYPE_1'))).filename;
%         % Name of the LTM  file containing the Ku-band CAL2 parameters
%         idx=find(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'CALIBRATION_TYPE_2',length('CALIBRATION_TYPE_2'))~=0);
%         if isempty(idx)
%            xref_altimeter_ltm_ku_cal2='Not available (not applied in the product)';
%         else
%            xref_altimeter_ltm_ku_cal2=files.sph.dsds(idx).filename;
%         end
%         % Name of the LTM file containing the C-band CAL2 parameters         
%         xref_altimeter_ltm_c_cal2='Not available for CR2';
%         % Name of the altimeter characterisation data file
%         xref_altimeter_characterisation=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'IPF_RA_DATABASE',length('IPF_RA_DATABASE'))).filename;
        % Semi-major axis of the reference ellipsoid
        semi_major_ellipsoid_axis=num2str(semi_major_axis_cst,15);
        % Flattening coeffcient of the reference ellipsoid
        ellipsoid_flattening=num2str(flat_coeff_cst,15);    
        % Attributes related to the OrbitalN_surfs_loc_estimatedormation required by Porto
%         orbit_phase_code=files.sph.orbit_info.phase_code;
%         orbit_cycle_num=((files.sph.orbit_info.cycle_num));
%         orbit_REL_Orbit=((files.sph.orbit_info.rel_orbit));
%         orbit_ABS_Orbit_Start=((files.sph.orbit_info.ABS_Orbit_Start));
%         orbit_Rel_Time_ASC_Node_Start=((files.sph.orbit_info.Rel_Time_ASC_Node_Start));
%         orbit_ABS_Orbit_Stop=((files.sph.orbit_info.ABS_Orbit_Stop));
%         orbit_Rel_Time_ASC_Node_Stop=((files.sph.orbit_info.Rel_Time_ASC_Node_Stop));
%         orbit_Equator_Cross_Time=(files.sph.orbit_info.Equator_Cross_Time);
%         orbit_Equator_Cross_Long=((files.sph.orbit_info.Equator_Cross_Long));
% 
end

%% PACKING L1B


%----------A. Time variables ----------------------------------------------
time_l1b_echo_sar_ku_name = 'time_l1b_echo_sar_ku';
id_aux      = netcdf.defVar(ncid,time_l1b_echo_sar_ku_name,double_type,ku_rec_dimension);
            netcdf.putAtt(ncid,id_aux,std_name_att,'time');
            netcdf.putAtt(ncid,id_aux,long_name_att,'UTC Seconds since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
            netcdf.putAtt(ncid,id_aux,calendar_name_att,'Gregorian');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
			netcdf.putAtt(ncid,id_aux,comment_att,'time at surface of the SAR measurement(multilooked waveform).');


UTC_day_l1b_echo_sar_ku_name = 'UTC_day_l1b_echo_sar_ku';
id_aux      = netcdf.defVar(ncid,UTC_day_l1b_echo_sar_ku_name,int16_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Days since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
            netcdf.putAtt(ncid,id_aux,units_att,day_units);
			netcdf.putAtt(ncid,id_aux,comment_att,'days elapsed since 2000-01-01. To be used to link with L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');


UTC_sec_l1b_echo_sar_ku_name = 'UTC_sec_l1b_echo_sar_ku';
id_aux 		= netcdf.defVar(ncid,UTC_sec_l1b_echo_sar_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,1.844674407370960e+19);
			netcdf.putAtt(ncid,id_aux,long_name_att,'Seconds in the day UTC, with microsecond resolution (Ku-band)');
			netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
			netcdf.putAtt(ncid,id_aux,comment_att,'seconds in the day. To be used to link with L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');



%----------B. Orbit and attitude variables ------------------------------
lat_l1b_echo_sar_ku_name = 'lat_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,lat_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude (positive N, negative S) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Latitude of measurement [-90, +90]: Positive at Nord, Negative at South');


lon_l1b_echo_sar_ku_name = 'lon_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,lon_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude (positive E, negative W) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'longitude of measurement [-180, +180]: Positive at East, Negative at West');


alt_l1b_echo_sar_ku_name = 'alt_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,alt_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite Centre of Mass');


orb_alt_rate_l1b_echo_sar_ku_name = 'orb_alt_rate_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,orb_alt_rate_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'orbital altitude rate');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Instantaneous altitude rate at the Centre of Mass');


%----------C. Flag time variables --------------------------------------

%----------D. Position/Velocity variables ------------------------------
x_pos_l1b_echo_sar_ku_name = 'x_pos_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,x_pos_l1b_echo_sar_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

y_pos_l1b_echo_sar_ku_name = 'y_pos_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,y_pos_l1b_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

z_pos_l1b_echo_sar_ku_name = 'z_pos_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,z_pos_l1b_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

x_vel_l1b_echo_sar_ku_name = 'x_vel_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,x_vel_l1b_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);

y_vel_l1b_echo_sar_ku_name = 'y_vel_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,y_vel_l1b_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);

z_vel_l1b_echo_sar_ku_name = 'z_vel_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,z_vel_l1b_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);

%----------D'. XYZ_Position of the measurement ------------------------------

x_surf_pos_l1b_echo_sar_ku_name = 'x_surf_pos_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,x_surf_pos_l1b_echo_sar_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Surface position-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

y_surf_pos_l1b_echo_sar_ku_name = 'y_surf_pos_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,y_surf_pos_l1b_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Surface position-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

z_surf_pos_l1b_echo_sar_ku_name = 'z_surf_pos_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,z_surf_pos_l1b_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Surface position-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);


satellite_mispointing_l1b_sar_echo_ku_name = 'satellite_mispointing_l1b_sar_echo_ku';
id_aux = netcdf.defVar(ncid,satellite_mispointing_l1b_sar_echo_ku_name,int32_type,[space_3D_dimension ,ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mispointing angle, measures by STRs: [1] Roll, [2] Pitch, [3] Yaw (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-7);
netcdf.putAtt(ncid,id_aux,comment_att,'Attitude mispointing, measured by STRs and post-processed by AOCS or by ground facility. The 3 components are given according to the ''space_3D'' dimension: [1] Roll, [2] Pitch, [3] Yaw. This variable includes the "mispointing bias" given by the variable mispointing_bias_ku. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: roll = pitch = yaw = 0');


%----------E. Navigation Bulletin ------------------------------

seq_count_l1b_echo_sar_ku_name = 'seq_count_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,seq_count_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,65535);
netcdf.putAtt(ncid,id_aux,long_name_att,'Sequence count');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%----------F. Operating instrument & tracking --------------------------
oper_instr_l1b_echo_sar_ku_name = 'oper_instr_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,oper_instr_l1b_echo_sar_ku_name,int8_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,127);
netcdf.putAtt(ncid,id_aux,long_name_att,'Operating instrument');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0b,1b');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'A,B (Sentinel-3) / Nominal, Redundant (CryoSat-2)');
netcdf.putAtt(ncid,id_aux,comment_att,'For Sentinel-3: Instrument A stands for SRAL Nominal and instrument B stands for SRAL Redundant');    

SAR_mode_l1b_echo_sar_ku_name='SAR_mode_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,SAR_mode_l1b_echo_sar_ku_name,int8_type, ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'SAR mode identifier');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0b,1b,2b (Sentinel-3) / 0b, 1b (Cryosat-2)');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'closed_loop, open_loop, open_loop_fixed_gain (Sentinel-3) / closed, open (CryoSat-2)');

%---------- G. H0, COR2 and AGC ----------------------------------------
h0_applied_l1b_echo_sar_ku_name = 'h0_applied_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,h0_applied_l1b_echo_sar_ku_name,uint32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,4294967295);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command H0');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d64d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'From ISP. Applied altitude command H0. Value the closest in time to the reference measurement');

cor2_applied_l1b_echo_sar_ku_name = 'cor2_applied_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,cor2_applied_l1b_echo_sar_ku_name,int16_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command COR2');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d1024d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'From ISP. Applied altitude variation. Value the closest in time to the reference measurement');


agccode_ku_l1b_echo_sar_ku_name = 'agccode_ku_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,agccode_ku_l1b_echo_sar_ku_name,int8_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'AGCCODE for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Computed: gain derived from the ATTCODE commands, using associated calibration values. Value the closest in time to the reference measurement');

%---------- H. Surface type -----------------------------------------------
surf_type_l1b_echo_sar_ku_name = 'surf_type_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,surf_type_l1b_echo_sar_ku_name,int8_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,127);
netcdf.putAtt(ncid,id_aux,long_name_att,'Altimeter surface type');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0,1,2,3');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'open_ocean or semi-enclosed_seas, enclosed_seas or lakes, continental_ice, land, Transponder');
netcdf.putAtt(ncid,id_aux,comment_att,'Surface Type at sub-satellite point. Value the closest in time to the reference measurement');

%---------- I. Altimeter range and Corrections ----------------------------
range_ku_l1b_echo_sar_ku_name = 'range_ku_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,range_ku_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Corrected range for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction');

uso_cor_l1b_echo_sar_ku_name = 'uso_cor_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,uso_cor_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'USO frequency drift correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

int_path_cor_ku_l1b_echo_sar_ku_name = 'int_path_cor_ku_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,int_path_cor_ku_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Internal path correction for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

if strcmp(mode,'SIN') && strcmp(mission,'CR2')&& strcmp(processing_mode_cnf,'SIN')
    int_path_2_cor_ku_l1b_echo_sar_ku_name = 'int_path_2_cor_ku_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,int_path_2_cor_ku_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Internal path correction for Ku band (second/receive only antenna CryoSat-2)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');
end

range_rate_l1b_echo_sar_ku_name = 'range_rate_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,range_rate_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Range rate');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%---------- J. AGC and Sigma0 scalings --------------------------------
scale_factor_ku_l1b_echo_sar_ku_name = 'scale_factor_ku_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,scale_factor_ku_l1b_echo_sar_ku_name,int32_type, [nl_dimension ku_rec_dimension]);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Scaling factor for sigma0 evaluation');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'This is a scaling factor in order to retrieve sigma-0 from the L1BS looks. It includes antenna gains and geometry satellite - surface. It is not applied to the L1BS looks');

%---------- K. Stack characterization--------------------------------------
% nb_stack_l1b_echo_sar_ku_name = 'nb_stack_l1b_echo_sar_ku';
% id_aux = netcdf.defVar(ncid,nb_stack_l1b_echo_sar_ku_name,uint16_type, ku_rec_dimension);
% %             netcdf.defVarFill(ncid,id_aux,false,65535);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Number of waveforms summed in stack');
% netcdf.putAtt(ncid,id_aux,units_att,number_units);
% 
look_angle_start_l1b_echo_sar_ku_name = 'look_angle_start_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,look_angle_start_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of first look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Look angle of the first contributing look (non-0 weight) to the L1B waveform');

look_angle_stop_l1b_echo_sar_ku_name = 'look_angle_stop_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,look_angle_stop_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of last look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Look angle of the last contributing look (non-0 weight) to the L1B waveform');

doppler_angle_start_l1b_echo_sar_ku_name = 'doppler_angle_start_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,doppler_angle_start_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of first look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Doppler angle of the first contributing look (non-0 weight) to the L1B waveform');

doppler_angle_stop_l1b_echo_sar_ku_name = 'doppler_angle_stop_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,doppler_angle_stop_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of last contributing look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Doppler angle of the last contributing look to the L1B waveform');

pointing_angle_start_l1b_echo_sar_ku_name = 'pointing_angle_start_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,pointing_angle_start_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of first contributing look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Pointing angle of the first contributing look to the L1B waveform');

pointing_angle_stop_l1b_echo_sar_ku_name = 'pointing_angle_stop_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,pointing_angle_stop_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of last contributing look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Pointing angle of the last contributing look to the L1B waveform');

burst_nb_l1bs_start_echo_sar_ku_name = 'burst_nb_l1bs_start_echo_sar_ku';
id_aux = netcdf.defVar(ncid,burst_nb_l1bs_start_echo_sar_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Burst index of first contributing look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Burst index of the first contributing look to the L1B waveform');

burst_nb_l1bs_stop_echo_sar_ku_name = 'burst_nb_l1bs_stop_echo_sar_ku';
id_aux = netcdf.defVar(ncid,burst_nb_l1bs_stop_echo_sar_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Burst index of last contributing look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Burst index of the last contributing look to the L1B waveform');


% 
% skew_stack_l1b_echo_sar_ku_name = 'skew_stack_l1b_echo_sar_ku';
% id_aux = netcdf.defVar(ncid,skew_stack_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
% %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Skewness of stack');
% netcdf.putAtt(ncid,id_aux,units_att,number_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
% netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'Skewness of the Gaussian that fits the integrated power of the looks within a stack. The skewness indicates how symmetric or asymmetric the power within the stack is');
% 
% kurt_stack_l1b_echo_sar_ku_name = 'kurt_stack_l1b_echo_sar_ku';
% id_aux = netcdf.defVar(ncid,kurt_stack_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
% %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Kurtosis of stack');
% netcdf.putAtt(ncid,id_aux,units_att,number_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
% netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'Kurtosis of the Gaussian that fits the integrated power of the looks within a stack. Kurtosis is a measure of peakiness');
% 
% stdev_stack_l1b_echo_sar_ku_name = 'stdev_stack_l1b_echo_sar_ku';
% id_aux = netcdf.defVar(ncid,stdev_stack_l1b_echo_sar_ku_name,uint32_type, ku_rec_dimension);
% %             netcdf.defVarFill(ncid,id_aux,false,4294967295);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: STD wrt look angle (Ku-band)');
% netcdf.putAtt(ncid,id_aux,units_att,rad_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
% netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'Standard deviation of the Gaussian that fits the integrated power of the looks within a stack. It is given with respect to the look angle. The width at -3dB of this Gaussian can be retrieved the following way: width_3db = 2*sqrt(2*ln2)*gaussian_fitting_std');
% 
% gaussian_fitting_centre_look_l1b_echo_sar_ku_name = 'gaussian_fitting_centre_look_l1b_echo_sar_ku';
% id_aux = netcdf.defVar(ncid,gaussian_fitting_centre_look_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: centre wrt look angle (Ku-band)');
% netcdf.putAtt(ncid,id_aux,units_att,rad_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
% netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'Position of the center of the Gaussian that fits the integrated power of the looks within a stack, with respect to look angle');
% 
% gaussian_fitting_centre_pointing_l1b_echo_sar_ku_name = 'gaussian_fitting_centre_pointing_l1b_echo_sar_ku';
% id_aux = netcdf.defVar(ncid,gaussian_fitting_centre_pointing_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: centre wrt pointing angle (Ku-band)');
% netcdf.putAtt(ncid,id_aux,units_att,rad_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
% netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'Position of the centre of the Gaussian that fits the integrated power of the looks within a stack, with respect to pointing angle');
% 
% beam_form_l1b_echo_sar_ku_name = 'beam_form_l1b_echo_sar_ku';
% id_aux = netcdf.defVar(ncid,beam_form_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Flag on beam formation quality in stack');
% netcdf.putAtt(ncid,id_aux,units_att,percent_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-2);
% netcdf.putAtt(ncid,id_aux,comment_att,'Position of the centre of the Gaussian that fits the integrated power of the looks within a stack, with respect to pointing angle');


%------------- L. Altimeter engineering variables -------------------------
altimeter_clock_l1b_echo_sar_ku_name = 'altimeter_clock_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,altimeter_clock_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Altimeter clock (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,Hz_units);
netcdf.putAtt(ncid,id_aux,add_offset_att,bw_ku_chd);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-9);
netcdf.putAtt(ncid,id_aux,comment_att,'This is the actual altimeter clock. The altimeter clock is based upon the USO clock. The nominal USO clock is 10MHz, while the nominal altimeter clock is 395 MHz. The actual USO clock is provided regularly as a driftN_surfs_loc_estimatedormation: based on the USO drift, the actual altimeter clock is computed, assuming a linear dependency');

pulse_repetition_interval_l1b_echo_sar_ku_name = 'pri_lrm_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,pulse_repetition_interval_l1b_echo_sar_ku_name,uint32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'PRI converted into seconds (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e19);
netcdf.putAtt(ncid,id_aux,comment_att,'The ''Pulse Repetition Interval''. PRI is constant within all received pulses in a radar cycle, but it can change within consecutive radar cycles. It is provided in counters of (T0) by the altimeter and converted in seconds by the L1 processor, with T0 = 1/altimeter_clock_ku');

%------------- M. Looks related variables -----------------------------
echoes_i_gprw_cor_l1bS_echo_sar_ku_name = 'echoes_i_gprw_cor_l1bS_echo_sar_ku';
id_aux = netcdf.defVar(ncid,echoes_i_gprw_cor_l1bS_echo_sar_ku_name,int8_type,[ns_dimension nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'I-samples for SAR L1B-S looks, arranged in stacks of NlxNs elements. I-samples are scaled to range [-127, +127] (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The i component of each L1B-S look. Each look is a fully calibrated, high resolution complex waveform. Each look within the stack is: (a) given in the frequency (range) domain, (b) aligned within the stack (slant range, Doppler range, window delay misalignments corrections applied), (c) fully calibrated. A final scaling, given in the variable ''q_scale_factor'', is applied in order to best fit the q component into 1 byte');

echoes_q_gprw_cor_l1bS_echo_sar_ku_name = 'echoes_q_gprw_cor_l1bS_echo_sar_ku';
id_aux = netcdf.defVar(ncid,echoes_q_gprw_cor_l1bS_echo_sar_ku_name,int8_type,[ns_dimension nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'Q-samples for SAR L1B-S looks, arranged in stacks of NlxNs elements. Q-samples are scaled to range [-127, +127] (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The q component of each L1B-S look. Each look is a fully calibrated, high resolution complex waveform. Each look within the stack is: (a) given in the frequency (range) domain, (b) aligned within the stack (slant range, Doppler range, window delay misalignments corrections applied), (c) fully calibrated. A final scaling, given in the variable ''q_scale_factor'', is applied in order to best fit the q component into 1 byte');

i_scale_factor_ku_name = 'i_scale_factor_ku_echo_sar_ku';
id_aux = netcdf.defVar(ncid,i_scale_factor_ku_name,float_type,[nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'I scale factor, to convert I samples from [-127, +127] to amplitude at antenna flange (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,sqrtW_per_count_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The i-samples scaling factor, computed in order to best fit the i-samples within 1 byte. The scaling, needed to convert the echoes_i_gprw_cor_l1bS_echo_sar_ku_name into sqrt(watt), is applied as follows: echoes_i_gprw_cor_l1bS_echo_sar_ku_sqr_watt(ku_rec,Nl, Ns) = echoes_i_gprw_cor_l1bS_echo_sar_ku(ku_rec,Nl, Ns) * i_scale_factor(ku_rec,Nl)');

q_scale_factor_ku_name = 'q_scale_factor_ku_echo_sar_ku';
id_aux = netcdf.defVar(ncid,q_scale_factor_ku_name,float_type,[nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'Q scale factor, to convert Q samples from [-127, +127] to amplitude at antenna flange (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,sqrtW_per_count_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The q-samples scaling factor, computed in order to best fit the q-samples within 1 byte. The scaling, needed to convert the echoes_q_gprw_cor_l1bS_echo_sar_ku_name into sqrt(watt), is applied as follows: echoes_q_gprw_cor_l1bS_echo_sar_ku_sqr_watt(ku_rec,Nl, Ns) = echoes_q_gprw_cor_l1bS_echo_sar_ku(ku_rec,Nl, Ns) * q_scale_factor(ku_rec,Nl)');

stack_mask_range_bin_l1b_echo_SAR_ku_name = 'stack_mask_range_bin_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,stack_mask_range_bin_l1b_echo_SAR_ku_name,uint8_type,[nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'Range bin stack mask ');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,zp_fact_range_cnf);
netcdf.putAtt(ncid,id_aux,comment_att,'Before the stack is multi-looked, the different looks are cropped according to these values, in order to clean the stack. For each look, the number of the first cropped sample is provided. The scale factor is equal to the range_oversampling_factor. The mask is not applied in the L1BS looks provided in this product.');

% waveform_scale_factor_l1b_echo_sar_ku_name = 'waveform_scale_factor_l1b_echo_sar_ku';
% id_aux = netcdf.defVar(ncid,waveform_scale_factor_l1b_echo_sar_ku_name,float_type,ku_rec_dimension);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Echo Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
% netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
% netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = power_waveform(ku_rec, Ns) * waveform_scale_factor(ku_rec)');
% 
zero_padding_l1b_echo_sar_ku_name = 'zero_padding_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,zero_padding_l1b_echo_sar_ku_name,int8_type,1);
netcdf.putAtt(ncid,id_aux,long_name_att,'Oversampling factor used in the range compression (FFT)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The instrument samples the waveforms with a 395 MHz clock, providing a nominal_sampling = c / 395e6 / 2 = ~0.379m (with c=speed of light). In addition, the ground processor can apply an oversampling factor, providing a waveform_sampling = nominal_sampling / range_oversampling_factor. Note that the altimeter range resolution is fixed and given by the chirp bandwidth of 320 MHz: c / 320e6 / 2 = ~0.468m');


%------------- N. Geophysical Corrections variables -----------------------

dry_tropo_correction_name = 'dry_tropo_correction_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,dry_tropo_correction_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Dry Tropospheric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

wet_tropo_correction_name = 'wet_tropo_correction_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,wet_tropo_correction_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Wet Tropospheric correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

inverse_baro_correction_name = 'inverse_baro_correction_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,inverse_baro_correction_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Inverse Barometric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

Dynamic_atmospheric_correction_name = 'Dynamic_atmospheric_correction_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,Dynamic_atmospheric_correction_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Dynamic Atmospheric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

GIM_iono_correction_name = 'GIM_iono_correction_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,GIM_iono_correction_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'GIM Ionospheric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

model_iono_correction_name = 'model_iono_correction_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,model_iono_correction_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Model Ionospheric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

ocean_equilibrium_tide_name = 'ocean_equilibrium_tide_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,ocean_equilibrium_tide_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Ocean Equilibrium Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

long_period_tide_height_name = 'long_period_tide_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,long_period_tide_height_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Long Period Ocean Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

ocean_loading_tide_name = 'ocean_loading_tide_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,ocean_loading_tide_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Ocean Loading Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

solid_earth_tide_name = 'solid_earth_tide_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,solid_earth_tide_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Solid Earth Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

geocentric_polar_tide_name = 'geocentric_polar_tide_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,geocentric_polar_tide_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Solid Earth Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%------------- O. Processing parameters used ------------------------------


%----------  Global Attributes definition -----------------------------------
%---- attributes inherited from Sentinel-3 product description-------------
id_aux = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,id_aux,'creation_time',date_creation);
netcdf.putAtt(ncid,id_aux,'Conventions',netcdf_v4_format);
netcdf.putAtt(ncid,id_aux,'mission_name',mission);
netcdf.putAtt(ncid,id_aux,'altimeter_sensor_name',altimeter_sensor_name);
netcdf.putAtt(ncid,id_aux,'gnss_sensor_name',gnss_sensor_name);
netcdf.putAtt(ncid,id_aux,'doris_sensor_name',doris_sensor_name);
netcdf.putAtt(ncid,id_aux,'doris_sensor_name',acq_station_name);
netcdf.putAtt(ncid,id_aux,'doris_sensor_name',acq_station_name);
netcdf.putAtt(ncid,id_aux,'first_meas_time',first_meas_time);
netcdf.putAtt(ncid,id_aux,'last_meas_time',last_meas_time);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_level0',xref_altimeter_level0);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_orbit',xref_altimeter_orbit);
netcdf.putAtt(ncid,id_aux,'xref_doris_USO',xref_doris_USO);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_sar_cal1',xref_altimeter_ltm_sar_cal1);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_ku_cal2',xref_altimeter_ltm_ku_cal2);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_c_cal2',xref_altimeter_ltm_c_cal2);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_characterisation',xref_altimeter_characterisation);
netcdf.putAtt(ncid,id_aux,'semi_major_ellipsoid_axis',semi_major_ellipsoid_axis);
netcdf.putAtt(ncid,id_aux,'ellipsoid_flattening',ellipsoid_flattening);
%--------------- add the attributes related to intermediate product--------
netcdf.putAtt(ncid,id_aux,'orbit_phase_code',orbit_phase_code);
netcdf.putAtt(ncid,id_aux,'orbit_cycle_num',orbit_cycle_num);
netcdf.putAtt(ncid,id_aux,'orbit_REL_Orbit',orbit_REL_Orbit);
netcdf.putAtt(ncid,id_aux,'orbit_ABS_Orbit_Start',orbit_ABS_Orbit_Start);
netcdf.putAtt(ncid,id_aux,'orbit_Rel_Time_ASC_Node_Start',orbit_Rel_Time_ASC_Node_Start);
netcdf.putAtt(ncid,id_aux,'orbit_ABS_Orbit_Stop',orbit_ABS_Orbit_Stop);
netcdf.putAtt(ncid,id_aux,'orbit_Rel_Time_ASC_Node_Stop',orbit_Rel_Time_ASC_Node_Stop);
netcdf.putAtt(ncid,id_aux,'orbit_Equator_Cross_Time',orbit_Equator_Cross_Time);
netcdf.putAtt(ncid,id_aux,'orbit_Equator_Cross_Long',orbit_Equator_Cross_Long);
netcdf.putAtt(ncid,id_aux,'orbit_Ascending_Flag',orbit_Ascending_Flag);

netcdf.endDef(ncid);

var_id = netcdf.inqVarID(ncid,'zero_padding_l1b_echo_sar_ku');
netcdf.putVar(ncid,var_id,0,zp_fact_range_cnf);

netcdf.close(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
