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
function [files] = create_NetCDF_L1A_S3(files, N_bursts)

global N_samples   semi_major_axis_cst flat_coeff_cst;
global bw_ku_chd zp_fact_range_cnf  sec_in_day_cst pi_cst
global mission mode N_max_beams_stack_chd processing_mode_cnf

% t6 = tic;

date_creation = datestr(now, '_yyyymmddTHHMMSS_');
switch mission
    case 'CR2'        
        switch mode
            case 'SAR'                
                files.filename_netCDF_BS = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA__A_',...
                                files.sph.product_info.product_id(20:20+30),...
                                date_creation,...
                                'isd','.nc');
            case 'SIN'                
                files.filename_netCDF_BS = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA__A_',...
                                files.sph.product_info.product_id(20:20+30),...
                                date_creation,...
                                'isd','.nc');
        end
    case {'S3','S3_','S3A','S3B'}
        files.filename_netCDF_BS = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA__A_',...
                                date_creation,...
                                'isd','.nc');
    case {'S6_'} 
end
ncid = netcdf.create('./results/data/measurement_l1a_reduced.nc','CLOBBER');

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

ku_rec_dimension = netcdf.defDim(ncid,'time_l1a_echo_sar_ku',N_bursts);
nl_dimension = netcdf.defDim(ncid,'Nl',64);
ns_dimension = netcdf.defDim(ncid,'Ns',N_samples);
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
rad_units = 'radian';
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
time_l1b_echo_sar_ku_name = 'time_l1a_echo_sar_ku';
id_aux      = netcdf.defVar(ncid,time_l1b_echo_sar_ku_name,double_type,ku_rec_dimension);
            netcdf.putAtt(ncid,id_aux,long_name_att,'UTC : l1a_echo_sar_ku mode');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);


seq_count_l1a_echo_sar_ku_name = 'seq_count_l1a_echo_sar_ku';
id_aux      = netcdf.defVar(ncid,seq_count_l1a_echo_sar_ku_name,int16_type, ku_rec_dimension);
            netcdf.putAtt(ncid,id_aux,'FillValue',65535);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Source sequence count : l1a_echo_sar_ku mode');
            netcdf.putAtt(ncid,id_aux,units_att,number_units);

UTC_day_l1a_echo_sar_ku_name = 'UTC_day_l1a_echo_sar_ku';
id_aux 		= netcdf.defVar(ncid,UTC_day_l1a_echo_sar_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,1.844674407370960e+19);
			netcdf.putAtt(ncid,id_aux,long_name_att,'day UTC : l1a_echo_sar_ku mode');
			netcdf.putAtt(ncid,id_aux,units_att,day_units);


UTC_sec_l1b_echo_sar_ku_name = 'UTC_sec_l1a_echo_sar_ku';
id_aux 		= netcdf.defVar(ncid,UTC_sec_l1b_echo_sar_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,1.844674407370960e+19);
			netcdf.putAtt(ncid,id_aux,long_name_att,'seconds in the day UTC : l1a_echo_sar_ku mode');
			netcdf.putAtt(ncid,id_aux,units_att,seconds_units);

uso_cor_l1b_echo_sar_ku_name = 'uso_cor_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,uso_cor_l1b_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'USO frequency drift correction : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

            
            
burst_count_prod_l1a_echo_sar_ku_name = 'burst_count_prod_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,burst_count_prod_l1a_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'bursts counter within the product : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'range 1 to number of records in the product');            
            

pitch_sat_pointing_l1a_echo_sar_ku_name = 'pitch_sat_pointing_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,pitch_sat_pointing_l1a_echo_sar_ku_name,int16_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'satellite pointing angle - pitch : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');

yaw_sat_pointing_l1a_echo_sar_ku_name = 'yaw_sat_pointing_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,yaw_sat_pointing_l1a_echo_sar_ku_name,int16_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'satellite pointing angle - yaw : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');

cog_cor_l1a_echo_sar_ku_name = 'cog_cor_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,cog_cor_l1a_echo_sar_ku_name,int16_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'Distance antenna-CoG correction : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Distance  in the z-component between the centre of mass of the satellite and the altimeter antenna reference point');



%----------D. Position/Velocity variables ------------------------------
x_pos_l1a_echo_sar_ku_name = 'x_pos_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,x_pos_l1a_echo_sar_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

y_pos_l1a_echo_sar_ku_name = 'y_pos_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,y_pos_l1a_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

z_pos_l1a_echo_sar_ku_name = 'z_pos_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,z_pos_l1a_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

x_vel_l1a_echo_sar_ku_name = 'x_vel_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,x_vel_l1a_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);

y_vel_l1a_echo_sar_ku_name = 'y_vel_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,y_vel_l1a_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);

z_vel_l1a_echo_sar_ku_name = 'z_vel_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,z_vel_l1a_echo_sar_ku_name,double_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);





alt_l1b_echo_sar_ku_name = 'alt_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,alt_l1b_echo_sar_ku_name,int32_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite Centre of Mass');


orb_alt_rate_l1b_echo_sar_ku_name = 'orb_alt_rate_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,orb_alt_rate_l1b_echo_sar_ku_name,int16_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'orbital altitude rate');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Instantaneous altitude rate at the Centre of Mass');

roll_sral_mispointing_l1a_echo_sar_ku_name = 'roll_sral_mispointing_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,roll_sral_mispointing_l1a_echo_sar_ku_name,int16_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'SRAL mispointing angle - roll : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');

pitch_sral_mispointing_l1a_echo_sar_ku_name = 'pitch_sral_mispointing_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,pitch_sral_mispointing_l1a_echo_sar_ku_name,int16_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'SRAL mispointing angle - pitch : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');

yaw_sral_mispointing_l1a_echo_sar_ku_name = 'yaw_sral_mispointing_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,yaw_sral_mispointing_l1a_echo_sar_ku_name,int16_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'SRAL mispointing angle - yaw : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');




%---------- G. H0, COR2 and AGC ----------------------------------------
h0_applied_l1a_echo_sar_ku_name = 'h0_applied_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,h0_applied_l1a_echo_sar_ku_name,uint32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,4294967295);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command H0');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d64d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'From ISP. Applied altitude command H0. Value the closest in time to the reference measurement');

cor2_applied_l1a_echo_sar_ku_name = 'cor2_applied_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,cor2_applied_l1a_echo_sar_ku_name,int16_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command COR2');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d1024d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'From ISP. Applied altitude variation. Value the closest in time to the reference measurement');

%---------- G. H0, COR2 and AGC ----------------------------------------
h0_nav_dem_l1a_echo_sar_ku_name = 'h0_nav_dem_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,h0_nav_dem_l1a_echo_sar_ku_name,uint32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,4294967295);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude command H0 computed with nav DEM : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d64d1e9_units);

cor2_nav_dem_l1a_echo_sar_ku_name = 'cor2_nav_dem_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,cor2_nav_dem_l1a_echo_sar_ku_name,int16_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude command COR2 computed with nav DEM : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d1024d1e9_units);


agc_ku_l1a_echo_sar_ku_name = 'agc_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,agc_ku_l1a_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected AGC for ku band: l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0.01);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'AGC corrected for instrumental errors (calibration)');


burst_count_cycle_l1a_echo_sar_ku_name = 'burst_count_cycle_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,burst_count_cycle_l1a_echo_sar_ku_name,int8_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'bursts counter within the cycle : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'range 1 to number of records in the product');            
           

sig0_cal_ku_l1a_echo_sar_ku_name = 'sig0_cal_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,sig0_cal_ku_l1a_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'internal calibration correction on Sigma0 for ku band: l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0.01);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);



%---------- I. Altimeter range and Corrections ----------------------------
range_ku_l1a_echo_sar_ku_name = 'range_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,range_ku_l1a_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Corrected range for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction');



int_path_cor_ku_l1a_echo_sar_ku_name = 'int_path_cor_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,int_path_cor_ku_l1a_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Internal path correction for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');



%---------- J. AGC and Sigma0 scalings --------------------------------
scale_factor_ku_l1a_echo_sar_ku_name = 'scale_factor_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,scale_factor_ku_l1a_echo_sar_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'scaling factor for sigma0 evaluation for ku band: l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'scaling factor corrected for AGC instrumental errors and internal calibration');


burst_power_cor_ku_l1a_echo_sar_ku_name = 'burst_power_cor_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,burst_power_cor_ku_l1a_echo_sar_ku_name,uint32_type,[nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'ku band burst power corrections (cal1) : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,fft_pow_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);

burst_phase_cor_ku_l1a_echo_sar_ku_name = 'burst_phase_cor_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,burst_phase_cor_ku_l1a_echo_sar_ku_name,uint32_type,[nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'ku band burst phase corrections (cal1) : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);

gprw_meas_ku_l1a_echo_sar_ku_name = 'gprw_meas_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,gprw_meas_ku_l1a_echo_sar_ku_name,uint32_type,[ns_dimension 3 ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'ku band samples of the normalized GPRW (cal2) : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,fft_pow_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);



%------------- M. Looks related variables -----------------------------
i_meas_ku_l1a_echo_sar_ku_name = 'i_meas_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,i_meas_ku_l1a_echo_sar_ku_name,int8_type,[ns_dimension nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'ku band echoes, i measurements : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'ku band echoes of a burst, I values (64*128 samples) in the time domain. The pulses are not corrected for AGC');

q_meas_ku_l1a_echo_sar_ku_name = 'q_meas_ku_l1a_echo_sar_ku';
id_aux = netcdf.defVar(ncid,q_meas_ku_l1a_echo_sar_ku_name,int8_type,[ns_dimension nl_dimension ku_rec_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'ku band echoes, q measurements : l1a_echo_sar_ku mode');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'ku band echoes of a burst, Q values (64*128 samples) in the time domain. The pulses are not corrected for AGC');


netcdf.endDef(ncid);


netcdf.close(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
