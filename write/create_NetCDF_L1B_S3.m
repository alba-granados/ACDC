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
% Last rev.: Monica Roca   / isardSAT 
% 
% Versions
% 1.0 
% 1.1 Updated time conversion for data between 2010 and 2016 (CR2)
% 2.0 Transformed to a function. Writting one record
% 2.1 zp changed to int16 as int8 did not work for zp > 128
% 2.2 case S3_ filename
% 2.3 case 'SIN' added 
% 2.4 added cnf flag processing_mode_cnf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function [files] = create_NetCDF_L1B_S3(files, N_surfs_loc_estimated)

global N_samples   semi_major_axis_cst flat_coeff_cst;
global bw_ku_chd zp_fact_range_cnf  sec_in_day_cst pi_cst
global mission mode N_max_beams_stack_chd
global compute_L1_stadistics_cnf include_wfms_aligned
global optional_ext_file_flag file_ext_string
%added by EM: 04.10.2016 
global ACDC_application_cnf cnf_p_ACDC
%global netcdf_type
global processing_mode_cnf
% t6 = tic;

date_creation = datestr(now, '_yyyymmddTHHMMSS_');

switch mission
    case 'CR2'        
        switch mode
            case 'SAR'
                files.filename_netCDF = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',...
                    files.sph.product_info.product_id(20:20+30),...
                    date_creation,...
                    'isd','.nc');
                if optional_ext_file_flag
                    files.filename_netCDF=strrep(files.filename_netCDF,'.nc',strcat(file_ext_string,'.nc'));
                end
            case 'SIN'                
                files.filename_netCDF = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',...
                                files.sph.product_info.product_id(20:20+30),...
                                date_creation,...
                                'isd','.nc');
        end
    case {'S3_','S3A','S3B'}
        files.filename_netCDF = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',...
                                date_creation,...
                                'isd','.nc');
%                 files.filename_netCDF = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',...
%                                 files.sph.product_info.product_id(20:20+30),...
%                                 date_creation,...
%                                 'isd','.nc');
    case {'S6_'} 
end
ncid = netcdf.create(files.filename_netCDF,'NETCDF4');
% switch netcdf_type
%     case 'netcdf4'
%         ncid = netcdf.create(files.filename_netCDF,'NETCDF4');
%     case 'netcdf3'
%         ncid = netcdf.create(files.filename_netCDF,'CLASSIC_MODEL');        
% end

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
% nl_dimension_size = N_max_beams_stack_chd;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;

ku_rec_dimension = netcdf.defDim(ncid,'time_l1b_echo_sar_ku',N_surfs_loc_estimated);
nl_dimension = netcdf.defDim(ncid,'max_multi_stack_ind',N_max_beams_stack_chd);
ns_dimension = netcdf.defDim(ncid,'echo_sample_ind',N_samples*zp_fact_range_cnf);
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
rad_units = 'rad';
percent_units='percent';


int8_type = 'NC_BYTE';
uint8_type = netcdf.getConstant('ubyte');
int16_type = 'NC_SHORT';
uint16_type = netcdf.getConstant('ushort');
int32_type = 'NC_INT';
uint32_type = netcdf.getConstant('uint');
int64_type = netcdf.getConstant('int64');
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
           [~,name_file,ext_file]=fileparts(files.characterization);
           xref_altimeter_ltm_ku_cal2=strcat(name_file,ext_file);
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
    case {'S3_','S3A','S3B'}
                %altimeter sensor name
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
%         orbit_Ascending_Flag=files.sph.orbit_info.Ascending_Flag;
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
id_aux = netcdf.defVar(ncid,UTC_sec_l1b_echo_sar_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,1.844674407370960e+19);
netcdf.putAtt(ncid,id_aux,long_name_att,'Seconds in the day UTC, with microsecond resolution (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,comment_att,'seconds in the day. To be used to link L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');


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


satellite_mispointing_l1b_sar_echo_ku_name = 'satellite_mispointing_l1b_sar_echo_ku';
id_aux = netcdf.defVar(ncid,satellite_mispointing_l1b_sar_echo_ku_name,int32_type,[space_3D_dimension ,ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mispointing angle, measures by STRs: [1] Roll, [2] Pitch, [3] Yaw (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-7);
netcdf.putAtt(ncid,id_aux,comment_att,'Attitude mispointing, measured by STRs and post-processed by AOCS or by ground facility. The 3 components are given according to the ''space_3D'' dimension: [1] Roll, [2] Pitch, [3] Yaw. This variable includes the "mispointing bias" given by the variable mispointing_bias_ku. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: roll = pitch = yaw = 0');

%----------C. Flag time variables --------------------------------------

%----------D. Position/Velocity variables ------------------------------
x_pos_l1b_echo_sar_ku_name = 'x_pos_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,x_pos_l1b_echo_sar_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-x component');
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

%----------E. Navigation Bulletin ------------------------------
seq_count_l1b_echo_sar_ku_name = 'seq_count_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,seq_count_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);        
% switch netcdf_type
%     case 'netcdf4'
%           id_aux = netcdf.defVar(ncid,seq_count_l1b_echo_sar_ku_name,uint16_type,ku_rec_dimension);        
%     case 'netcdf3'
%            %increase the size
%           id_aux = netcdf.defVar(ncid,seq_count_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
% end    
%             netcdf.defVarFill(ncid,id_aux,false,65535);
netcdf.putAtt(ncid,id_aux,long_name_att,'Sequence count');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%----------F. Operating instrument & tracking --------------------------
oper_instr_l1b_echo_sar_ku_name = 'oper_instr_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,oper_instr_l1b_echo_sar_ku_name,int8_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,127);
netcdf.putAtt(ncid,id_aux,long_name_att,'Operating instrument');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0b,1b');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'A,B (Sentinel-3) / Nominal, Redundant (CryoSat-2)');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement. For Sentinel-3: Instrument A stands for SRAL Nominal and instrument B stands for SRAL Redundant');    

SAR_mode_l1b_echo_sar_ku_name='SAR_mode_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,SAR_mode_l1b_echo_sar_ku_name,int8_type, ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'SAR mode identifier');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0b,1b,2b (Sentinel-3) / 0b, 1b (Cryosat-2)');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'closed_loop, open_loop, open_loop_fixed_gain (Sentinel-3) / closed, open (CryoSat-2)');

%---------- G. H0, COR2 and AGC ----------------------------------------
h0_applied_l1b_echo_sar_ku_name = 'h0_applied_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,h0_applied_l1b_echo_sar_ku_name,int64_type, ku_rec_dimension);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,h0_applied_l1b_echo_sar_ku_name,int64_type, ku_rec_dimension);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,h0_applied_l1b_echo_sar_ku_name,int64_type, ku_rec_dimension);
% end
%             netcdf.defVarFill(ncid,id_aux,false,4294967295);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command H0');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d64d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

cor2_applied_l1b_echo_sar_ku_name = 'cor2_applied_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,cor2_applied_l1b_echo_sar_ku_name,int16_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command COR2');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d1024d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');


agccode_ku_l1b_echo_sar_ku_name = 'agccode_ku_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,agccode_ku_l1b_echo_sar_ku_name,int8_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'AGCCODE for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%---------- H. Surface type -----------------------------------------------
surf_type_l1b_echo_sar_ku_name = 'surf_type_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,surf_type_l1b_echo_sar_ku_name,int8_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,127);
netcdf.putAtt(ncid,id_aux,long_name_att,'Altimeter surface type');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0,1,2,3');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'open_ocean or semi-enclosed_seas, enclosed_seas or lakes, continental_ice, land,Transponder');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%---------- I. Altimeter range and Corrections ----------------------------
range_ku_l1b_echo_sar_ku_name = 'range_ku_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,range_ku_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Corrected range for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction');

if include_wfms_aligned
%EM 14.04.2016: window delay alignment
range_ku_surf_aligned_l1b_echo_sar_ku_name = 'range_ku_surf_aligned_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,range_ku_surf_aligned_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Corrected range for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction, including the additional correction for leading edge alginment w.r.t first surface');
end

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
id_aux = netcdf.defVar(ncid,scale_factor_ku_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Scaling factor for sigma0 evaluation');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'This is a scaling factor in order to retrieve sigma-0 from the L1B waveform. It includes antenna gains and geometry satellite - surface. It is not applied to the L1B waveforms');

%---------- K. Stack characterization--------------------------------------
nb_stack_l1b_echo_sar_ku_name = 'nb_stack_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,nb_stack_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,nb_stack_l1b_echo_sar_ku_name,uint16_type, ku_rec_dimension);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,nb_stack_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
% end
%             netcdf.defVarFill(ncid,id_aux,false,65535);
netcdf.putAtt(ncid,id_aux,long_name_att,'Number of waveforms summed in stack (contributing beams: effective number of looks or beams from stack used in multilooking; if beams with all samples set to zero not to be included in multilooking they are accordingly not accounted in this number; if there are some gaps in-between mask set to zero all the samples related to those beams and they will not be contributing to the multilooking )');
netcdf.putAtt(ncid,id_aux,units_att,number_units);

nb_stack_start_stop_l1b_echo_sar_ku_name = 'nb_stack_start_stop_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,nb_stack_start_stop_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Number of waveforms in stack (considering the number of beams/looks from start and stop beams: they correspond to the first and last beams at edges of stack. If option of discarding beams with all-zeros samples these correspond to the first and last beams with non-all zeros samples (if there exist gaps in between: mask setting to zero all samples, these in-between beams are considered anyway in the "nb_stack_start_stop_l1b_echo_sar_ku"), otherwise they correspond to the very first and very last beam in the stack ). This number of beams will be useful to construct accordingly the modelled stack for retracking.');
netcdf.putAtt(ncid,id_aux,units_att,number_units);

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

if compute_L1_stadistics_cnf
    skew_stack_l1b_echo_sar_ku_name = 'skew_stack_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,skew_stack_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Skewness of stack');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Skewness of the Gaussian that fits the integrated power of the looks within a stack. The skewness indicates how symmetric or asymmetric the power within the stack is');
    
    kurt_stack_l1b_echo_sar_ku_name = 'kurt_stack_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,kurt_stack_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Kurtosis of stack');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Kurtosis of the Gaussian that fits the integrated power of the looks within a stack. Kurtosis is a measure of peakiness');
    
    stdev_stack_l1b_echo_sar_ku_name = 'stdev_stack_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,stdev_stack_l1b_echo_sar_ku_name,int64_type, ku_rec_dimension);
    % switch netcdf_type
    %     case 'netcdf4'
    %         id_aux = netcdf.defVar(ncid,stdev_stack_l1b_echo_sar_ku_name,uint32_type, ku_rec_dimension);
    %     case 'netcdf3'
    %         id_aux = netcdf.defVar(ncid,stdev_stack_l1b_echo_sar_ku_name,int64_type, ku_rec_dimension);
    % end
    %             netcdf.defVarFill(ncid,id_aux,false,4294967295);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: STD wrt look angle (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Standard deviation of the Gaussian that fits the integrated power of the looks within a stack. It is given with respect to the look angle. The width at -3dB of this Gaussian can be retrieved the following way: width_3db = 2*sqrt(2*ln2)*gaussian_fitting_std');
    
    gaussian_fitting_centre_look_l1b_echo_sar_ku_name = 'gaussian_fitting_centre_look_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,gaussian_fitting_centre_look_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: centre wrt look angle (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Position of the center of the Gaussian that fits the integrated power of the looks within a stack, with respect to look angle');
    
    gaussian_fitting_centre_pointing_l1b_echo_sar_ku_name = 'gaussian_fitting_centre_pointing_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,gaussian_fitting_centre_pointing_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: centre wrt pointing angle (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Position of the centre of the Gaussian that fits the integrated power of the looks within a stack, with respect to pointing angle');
    
    beam_form_l1b_echo_sar_ku_name = 'beam_form_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,beam_form_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Flag on beam formation quality in stack');
    netcdf.putAtt(ncid,id_aux,units_att,percent_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-2);
    netcdf.putAtt(ncid,id_aux,comment_att,'Beam formation quality in percentage');
end


%------------- L. Altimeter engineering variables -------------------------
altimeter_clock_l1b_echo_sar_ku_name = 'altimeter_clock_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,altimeter_clock_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Altimeter clock (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,Hz_units);
netcdf.putAtt(ncid,id_aux,add_offset_att,bw_ku_chd);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-9);
netcdf.putAtt(ncid,id_aux,comment_att,'This is the actual altimeter clock.');

pulse_repetition_interval_l1b_echo_sar_ku_name = 'pri_lrm_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,pulse_repetition_interval_l1b_echo_sar_ku_name,int64_type,ku_rec_dimension);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,pulse_repetition_interval_l1b_echo_sar_ku_name,uint32_type,ku_rec_dimension);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,pulse_repetition_interval_l1b_echo_sar_ku_name,int64_type,ku_rec_dimension);
% end
netcdf.putAtt(ncid,id_aux,long_name_att,'PRI converted into seconds (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-12);
netcdf.putAtt(ncid,id_aux,comment_att,'The ''Pulse Repetition Interval''. PRI is constant within all received pulses in a radar cycle, but it can change within consecutive radar cycles.');

%------------- M. Waveform related variables -----------------------------
i2q2_meas_ku_l1b_echo_sar_ku_name = 'i2q2_meas_ku_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,i2q2_meas_ku_l1b_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,i2q2_meas_ku_l1b_echo_sar_ku_name,uint16_type,[ns_dimension ku_rec_dimension]);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,i2q2_meas_ku_l1b_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
% end
netcdf.putAtt(ncid,id_aux,long_name_att,'SAR	Power Echo waveform: scaled	0-65535 (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The SAR L1B Power waveforms is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo_sar_ku'', in order to best fit the waveform into 2 bytes');

if include_wfms_aligned
%EM: 14.04.2016
i2q2_meas_ku_wdcorr_l1b_echo_sar_ku_name = 'i2q2_meas_ku_wdcorr_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,i2q2_meas_ku_wdcorr_l1b_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,i2q2_meas_ku_wdcorr_l1b_echo_sar_ku_name,uint16_type,[ns_dimension ku_rec_dimension]);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,i2q2_meas_ku_wdcorr_l1b_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
% end
netcdf.putAtt(ncid,id_aux,long_name_att,'SAR	Power Echo waveform aligned surfaces: scaled	0-65535 (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The SAR L1B Power waveforms (aligned w.r.t reference surface) is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo_sar_ku'', in order to best fit the waveform into 2 bytes. The waveforms have been algined taking into account window delay correction w.r.t first surface.');
end

stack_mask_range_bin_l1b_echo_sar_ku_name = 'stack_mask_range_bin_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,stack_mask_range_bin_l1b_echo_sar_ku_name,uint8_type,[nl_dimension ku_rec_dimension]);
%netcdf.defVarFill(ncid,id_aux,false,0);
netcdf.putAtt(ncid,id_aux,long_name_att,'Range bin stack mask ');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,zp_fact_range_cnf);
netcdf.putAtt(ncid,id_aux,comment_att,'The zero-mask applied to the stack before multilooking. Each element of the mask refers to a look in the stack and indicates the index of the first sample set to zero. The first num_looks_start_stop elements of the mask are valid, while the remaining ones are filled with 0.');

waveform_scale_factor_l1b_echo_sar_ku_name = 'waveform_scale_factor_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,waveform_scale_factor_l1b_echo_sar_ku_name,float_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Echo Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_l1b_echo_sar_ku(ku_rec, Ns) * waveform_scale_factor_l1b_echo_sar_ku(ku_rec)');

zero_padding_l1b_echo_sar_ku_name = 'zero_padding_l1b_echo_sar_ku';
id_aux = netcdf.defVar(ncid,zero_padding_l1b_echo_sar_ku_name,int16_type,[]); %scalar
netcdf.putAtt(ncid,id_aux,long_name_att,'Oversampling factor used in the range compression (FFT)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'the ground processor can apply an oversampling factor, providing a waveform_sampling = nominal_sampling / range_oversampling_factor. Note that the altimeter range resolution is fixed and given by the chirp bandwidth');

%EM: 22.09.2016
if strcmp(mode,'SIN') && strcmp(mission,'CR2') && strcmp(processing_mode_cnf,'SIN') 
    i2q2_meas_ku_l1b_echo_sar_ku_2_name = 'i2q2_meas_ku_l1b_echo_sar_ku_2';
    id_aux = netcdf.defVar(ncid,i2q2_meas_ku_l1b_echo_sar_ku_2_name,int32_type,[ns_dimension ku_rec_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'SAR Power Echo waveform 2: scaled	0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The SAR L1B Power waveform for Channel 2 (SARin). It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo_sar_ku'', in order to best fit the waveform into 2 bytes');
    
    waveform_scale_factor_l1b_echo_sar_ku_2_name = 'waveform_scale_factor_l1b_echo_sar_ku_2';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_l1b_echo_sar_ku_2_name,float_type,ku_rec_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Echo Scale Factor Channel 2, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform for Channel 2 (SARin) scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_l1b_echo_sar_ku(ku_rec, Ns) * waveform_scale_factor_l1b_echo_sar_ku(ku_rec)');
    
    coherence_meas_ku_l1b_echo_sar_ku_name = 'coherence_meas_ku_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,coherence_meas_ku_l1b_echo_sar_ku_name,int16_type,[ns_dimension ku_rec_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Interferometric Coherence');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Coherence between the two interferometric channels operating in SARIn, both channels are fully calibrated.');
    
    phase_difference_meas_ku_l1b_echo_sar_ku_name = 'phase_difference_meas_ku_l1b_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,phase_difference_meas_ku_l1b_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Phase difference SARIn  (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,comment_att,'Interferometric Phase multilooked obtained from the pair of channels operating in SARIn');
   
    
    if include_wfms_aligned
        
        i2q2_meas_ku_wdcorr_l1b_echo_sar_ku_2_name = 'i2q2_meas_ku_wdcorr_l1b_echo_sar_ku_2';
        id_aux = netcdf.defVar(ncid,i2q2_meas_ku_wdcorr_l1b_echo_sar_ku_2_name,int32_type,[ns_dimension ku_rec_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'SAR Power Echo waveform channel 2 aligned surfaces: scaled	0-65535 (Ku-band)');
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'The SAR L1B Power waveform channel 2 (aligned w.r.t given reference surface). It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo_sar_ku'', in order to best fit the waveform into 2 bytes. The waveforms have been algined taking into account window delay correction w.r.t first surface.');
        
        coherence_meas_ku_wdcorr_l1b_echo_sar_ku_name = 'coherence_meas_ku_wdcorr_l1b_echo_sar_ku';
        id_aux = netcdf.defVar(ncid,coherence_meas_ku_wdcorr_l1b_echo_sar_ku_name,int16_type,[ns_dimension ku_rec_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'Interferometric Coherence (aligned surfaces)');
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
        netcdf.putAtt(ncid,id_aux,comment_att,'Coherence between the two interferometric channels operating in SARIn (aligned w.r.t given reference surface), both channels are fully calibrated.');
        
        phase_difference_meas_ku_wdcorr_l1b_echo_sar_ku_name = 'phase_difference_meas_ku_wdcorr_l1b_echo_sar_ku';
        id_aux = netcdf.defVar(ncid,phase_difference_meas_ku_wdcorr_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
        netcdf.putAtt(ncid,id_aux,long_name_att,'Phase difference SARIn aligned surfaces  (Ku-band)');
        netcdf.putAtt(ncid,id_aux,units_att,rad_units);
        netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
        netcdf.putAtt(ncid,id_aux,comment_att,'Interferometric Phase multilooked obtained from the pair of channels operating in SARIn (aligned w.r.t given reference surface)');

    end
           
end

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
netcdf.putAtt(ncid,id_aux,long_name_att,'Geocentric Polar Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%------------- O. Processing parameters used ------------------------------


%added EM 04.10.2016
% %------------- P. ACDC variables ----------------------------------------
if ACDC_application_cnf
    %-------------- Waveform -----------------------------------------------
    %multilooked ACDC waveform
    i2q2_meas_ku_ACDC_echo_sar_ku_name = 'i2q2_meas_ku_ACDC_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,i2q2_meas_ku_ACDC_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC SAR Power Echo waveform: scaled 0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The ACDC SAR L1B Power waveform. It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_ACDC_echo_sar_ku'', in order to best fit the waveform into 2 bytes');
    %fitted ACDC waveform
    i2q2_fit_ku_ACDC_echo_sar_ku_name = 'i2q2_fit_ku_ACDC_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,i2q2_fit_ku_ACDC_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC SAR Fitted Power Echo waveform: scaled 0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The fitted ACDC SAR L1B Power waveform. It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_ACDC_echo_sar_ku'', in order to best fit the waveform into 2 bytes');
    %fitted L1B SAR waveform
    i2q2_fit_ku_echo_sar_ku_name = 'i2q2_fit_ku_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,i2q2_fit_ku_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'SAR Fitted Power Echo waveform: scaled 0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The fitted SAR L1B Power waveform. It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_ACDC_echo_sar_ku'', in order to best fit the waveform into 2 bytes');
        
    
    %-------------- Scaling -----------------------------------------------
    %scaling ACDC waveform
    waveform_scale_factor_ACDC_echo_sar_ku_name = 'waveform_scale_factor_ACDC_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_ACDC_echo_sar_ku_name,float_type,ku_rec_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Echo Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The ACDC waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_ACDC_echo_sar_ku(ku_rec, Ns) * waveform_scale_factor_ACDC_echo_sar_ku(ku_rec)');
    %scaling fitted ACDC waveform
    waveform_scale_factor_fit_ACDC_echo_sar_ku_name = 'waveform_scale_factor_fit_ACDC_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_fit_ACDC_echo_sar_ku_name,float_type,ku_rec_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC fit Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The ACDC fitted waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_ACDC_echo_sar_ku(ku_rec, Ns) * waveform_scale_factor_ACDC_echo_sar_ku(ku_rec)');    
    %scaling fitted ACDC waveform
    waveform_scale_factor_fit_echo_sar_ku_name = 'waveform_scale_factor_fit_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_fit_echo_sar_ku_name,float_type,ku_rec_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'SAR fit Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The SAR fitted waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_ACDC_echo_sar_ku(ku_rec, Ns) * waveform_scale_factor_ACDC_echo_sar_ku(ku_rec)');
    %-------------- Range index -------------------------------------------
    range_index_ACDC_echo_sar_ku_name = 'range_index_ACDC_echo_sar_ku';
    id_aux = netcdf.defVar(ncid,range_index_ACDC_echo_sar_ku_name,int32_type,[ns_dimension ku_rec_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC range index (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'The ACDC range index for ACDC SAR Power waveform.');
    
    %-------------- Geophysical retrievals --------------------------------
    %******************* retracked range **********************************
    retracked_range_ACDC_20_ku_name = 'retracked_range_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,retracked_range_ACDC_20_ku_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Retracked range for Ku band (2-step analytical retracker)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC corrected range by the retracker offset, the reference range includes instrumental corrections (already the USO frequency drift and the internal/instrument corrections).');
    %******************* epoch ********************************************
    epoch_ACDC_20_ku_name = 'epoch_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,epoch_ACDC_20_ku_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Epoch (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Estimated ACDC epoch in seconds w.r.t center of the window (window delay is given to the center of the window) using the 2-step analytical retracker. This corresponds to zero-padded sample value.');
    %******************** SWH  ********************************************
    swh_ACDC_20_ku_name = 'swh_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,swh_ACDC_20_ku_name,int16_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Significant waveheight (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Fitted significant waveheight. This corresponds to 4 times the fitted standard deviation of the surface height.');
    swh_ACDC_PRE_20_ku_name = 'swh_ACDC_PRE_20_ku';
    id_aux = netcdf.defVar(ncid,swh_ACDC_PRE_20_ku_name,int16_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Pre-retracker ACDC Significant waveheight (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Pre-retracker ACDC Fitted significant waveheight. This corresponds to 4 times the fitted standard deviation of the surface height.');
    %******************** SSH  ********************************************
    ssh_ACDC_20_ku_name = 'ssh_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,ssh_ACDC_20_ku_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Surface height (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
    ssh_ACDC_PRE_20_ku_name = 'ssh_ACDC_PRE_20_ku';
    id_aux = netcdf.defVar(ncid,ssh_ACDC_PRE_20_ku_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Pre-retracker ACDC Surface height (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Pre-retracker ACDC Surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
    %******************** SIGMA0 ******************************************    
    sig0_ACDC_20_ku_name = 'sig0_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,sig0_ACDC_20_ku_name,int16_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,32767);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Backscattering coefficient (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
    
    sig0_ACDC_PRE_20_ku_name = 'sig0_ACDC_PRE_20_ku';
    id_aux = netcdf.defVar(ncid,sig0_ACDC_PRE_20_ku_name,int16_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,32767);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Pre-retracker ACDC Backscattering coefficient (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Pre-retracker ACDC Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
    
    %********************* Pu: fitted peak power **************************    
    Pu_analytical_20_ku_name = 'Pu_analytical_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,Pu_analytical_20_ku_name,int16_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Fitted peak power (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Peak power of the fitted waveform.');
    %********************** Pearson Correlation Coefficient ***************
    Pearson_corr_ACDC_20_ku_name = 'Pearson_corr_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,Pearson_corr_ACDC_20_ku_name,int16_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Pearson coefficient (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,percent_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Pearson correlation coefficient as percentage indicating the goodness of fitting between the real ACDC waveform and the fitted one.');

    Pearson_corr_ACDC_PRE_20_ku_name = 'Pearson_corr_ACDC_PRE_20_ku';
    id_aux = netcdf.defVar(ncid,Pearson_corr_ACDC_PRE_20_ku_name,int16_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Pre-retracker ACDC Pearson coefficient (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,percent_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Pre-retracker ACDC Pearson correlation coefficient as percentage indicating the goodness of fitting between the real ACDC waveform and the fitted one.');    
    
    %************************ Fitting flag ********************************    
    Flag_fitting_ACDC_20_ku_name = 'Flag_fitting_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,Flag_fitting_ACDC_20_ku_name,int8_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,127);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Exit flag fitting');
    netcdf.putAtt(ncid,id_aux,flag_values_att,'1,2,3,4,0,-1,-2');
    netcdf.putAtt(ncid,id_aux,flag_desc_att,'1: Function converged to a solution x; 2: Change in x was less than the specified tolerance; 3: Change in the residual was less than the specified tolerance; 4: Magnitude of search direction was smaller than the specified tolerance; 0: Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations; -1: Output function terminated the algorithm; -2: Problem is infeasible: the bounds lb and ub are inconsistent');
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Flag on the reason the LSE solver stopped.');
      
end


%----------  Global Attributes definition -----------------------------------
%---- attributes inherited from Sentinel-3 product description-------------
id_aux = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,id_aux,'creation_time',date_creation);
netcdf.putAtt(ncid,id_aux,'Conventions',netcdf_v4_format);
netcdf.putAtt(ncid,id_aux,'mission_name',mission);
netcdf.putAtt(ncid,id_aux,'altimeter_sensor_name',altimeter_sensor_name);
netcdf.putAtt(ncid,id_aux,'gnss_sensor_name',gnss_sensor_name);
netcdf.putAtt(ncid,id_aux,'doris_sensor_name',doris_sensor_name);
netcdf.putAtt(ncid,id_aux,'acq_station_name',acq_station_name);
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
netcdf.putVar(ncid,var_id,zp_fact_range_cnf);

netcdf.close(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
