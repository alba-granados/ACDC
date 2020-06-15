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
% 2.1 zp changed to int16 as int8 did not work for zp > 128
% 2.2 case S3_ filename
% 2.3 case 'SIN' added
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function L2 = create_NetCDF_L2_S3 (L2,SPH,mission,mode,output_path_L2_ISD)


% input_path_L1B_ISD = [main_path,'L1B\isardSAT\'];
% input_path_L1B_ESA = [main_path,'L1B\ESA\'];
% input_path_L2_ESA = [main_path,'L2\ESA\'];
% path_comparison_results = [main_path,'L2\comparison\'];

global semi_major_axis_cst flat_coeff_cst

%global netcdf_type
% t6 = tic;

date_creation = datestr(now, '_yyyymmddTHHMMSS_');
switch mission
    case 'CR2'        
        switch mode
            case 'SAR'                
                L2.filename_netCDF = strcat(strcat(output_path_L2_ISD),mission,'_SR_1_SRA____',...
                                SPH.product_info.product_id(20:20+30),...
                                date_creation,...
                                'isd','.nc');
            case 'SIN'                
                L2.filename_netCDF = strcat(strcat(output_path_L2_ISD),mission,'_SR_1_SRA____',...
                                SPH.product_info.product_id(20:20+30),...
                                date_creation,...
                                'isd','.nc');
        end
    case {'S3_','S3A','S3B'}
        L2.filename_netCDF = strcat(strcat(output_path_L2_ISD,'data/'),mission,'_SR_1_SRA____',...
                                date_creation,...
                                'isd','.nc');
%                 files.filename_netCDF = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',...
%                                 files.sph.product_info.product_id(20:20+30),...
%                                 date_creation,...
%                                 'isd','.nc');
    case {'S6_'} 
end

L2.ncid = netcdf.create(L2.filename_netCDF,'NETCDF4');

% switch netcdf_type
%     case 'netcdf4'
%         L2.ncid = netcdf.create(files.filename_netCDF,'NETCDF4');
%     case 'netcdf3'
%         L2.ncid = netcdf.create(files.filename_netCDF,'CLASSIC_MODEL');        
% end

long_name_att = 'long_name';
std_name_att = 'standard_name';
calendar_name_att='calendar';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';

netcdf_v4_format = 'netcdf4';

% ku_rec_dimension = 'time_l1b_echo_sar_ku';
% nl_dimension = 'Nl';
% nl_dimension_size = N_max_beams_stack_chd;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;

ku_rec_dimension = netcdf.defDim(L2.ncid,'time_20Hz_ku',L2.ISD_num_surfaces_filtered);
single_dimension = netcdf.defDim(L2.ncid,'single',1);

day_units = 'day';
seconds_units = 'seconds';

degrees_units = 'degrees';
meters_units = 'meters';
dB_units = 'dB';

int16_type = 'NC_SHORT';
int32_type = 'NC_INT';
double_type= 'NC_DOUBLE';


%% CODING L2

%--------- Global Attribtues of the netCDF ---------------------------
switch mission
    case 'CR2'
        %altimeter sensor name
        switch SPH.ins_info.ins_id
            case 'A'
               altimeter_sensor_name='SIRAL Nominal'; 
            case 'B'
               altimeter_sensor_name='SIRAL Redundant'; 
        end
        altimeter_sensor_name=deblank(SPH.ins_info.ins_conf);
        
        % UTC date of the first measurement
        first_meas_time=SPH.product_info.product_time_info.produc_start_time;
        
        % UTC date of the last measurement
        last_meas_time=SPH.product_info.product_time_info.produc_stop_time;
        
        % Value of the first valid latitude
        first_meas_lat=L2.ISD_lat_surf(1)*1e-6;
        % Value of the last valid latitude
        last_meas_lat=L2.ISD_lat_surf(end)*1e-6;
        % Value of the first valid longitude
        first_meas_lon=L2.ISD_lon_surf(1)*1e-6;
        % Value of the last valid longitude
        last_meas_lon=L2.ISD_lon_surf(end)*1e-6;
        
        % Semi-major axis of the reference ellipsoid
        semi_major_ellipsoid_axis=num2str(semi_major_axis_cst,15);
        
        % Flattening coeffcient of the reference ellipsoid
        ellipsoid_flattening=num2str(flat_coeff_cst,15);    
        
        % Attributes related to the Orbital information required by Porto
        orbit_cycle_num=((SPH.orbit_info.cycle_num));
        orbit_REL_Orbit=((SPH.orbit_info.rel_orbit));
        orbit_ABS_Orbit_Start=((SPH.orbit_info.ABS_Orbit_Start));
        
    
    
    case {'S3_','S3A','S3B'}
        %altimeter sensor name
        altimeter_sensor_name='SRAL '; 
        
        % UTC date of the first measurement
        first_meas_time=SPH.product_info.product_time_info.produc_start_time;
        % UTC date of the last measurement
        last_meas_time=SPH.product_info.product_time_info.produc_stop_time;
        
        % Value of the first valid latitude
        first_meas_lat=L2.ISD_lat_surf(1);
        % Value of the last valid latitude
        last_meas_lat=L2.ISD_lat_surf(last);
        % Value of the first valid longitude
        first_meas_lon=L2.ISD_lon_surf(1);
        % Value of the last valid longitude
        last_meas_lon=L2.ISD_lon_surf(last);
        
        % Semi-major axis of the reference ellipsoid
        semi_major_ellipsoid_axis=num2str(semi_major_axis_cst,15);
        
        % Flattening coeffcient of the reference ellipsoid
        ellipsoid_flattening=num2str(flat_coeff_cst,15);    
        

        % Attributes related to the Orbital information required by Porto
        orbit_cycle_num=((SPH.orbit_info.cycle_num));
        orbit_REL_Orbit=((SPH.orbit_info.rel_orbit));
        orbit_ABS_Orbit_Start=((SPH.orbit_info.ABS_Orbit_Start));
        
end

%% PACKING L2


%---------- Time variables ----------------------------------------------
time_20Hz_ku_name = 'time_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,time_20Hz_ku_name,double_type,ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,std_name_att,'time');
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'UTC Seconds since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
netcdf.putAtt(L2.ncid,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(L2.ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'time at surface of the SAR measurement(multi-looked waveform).');

UTC_day_20Hz_ku_name = 'UTC_day_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,UTC_day_20Hz_ku_name,int16_type, ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Days since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
netcdf.putAtt(L2.ncid,id_aux,units_att,day_units);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'days elapsed since 2000-01-01. To be used to link with L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');

UTC_sec_20Hz_ku_name = 'UTC_sec_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,UTC_sec_20Hz_ku_name,double_type,ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Seconds in the day UTC, with microsecond resolution (Ku-band)');
netcdf.putAtt(L2.ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'seconds in the day. To be used to link L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');



%---------- Orbit variables ------------------------------
lat_20Hz_ku_name = 'lat_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,lat_20Hz_ku_name,int32_type,ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'latitude (positive N, negative S) (Ku-band)');
netcdf.putAtt(L2.ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'Latitude of measurement [-90, +90]: Positive at Nord, Negative at South');


lon_20Hz_ku_name = 'lon_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,lon_20Hz_ku_name,int32_type,ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'longitude (positive E, negative W) (Ku-band)');
netcdf.putAtt(L2.ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'longitude of measurement [-180, +180]: Positive at East, Negative at West');


alt_20Hz_ku_name = 'alt_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,alt_20Hz_ku_name,int32_type,ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'altitude of satellite');
netcdf.putAtt(L2.ncid,id_aux,units_att,meters_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'Altitude of the satellite Centre of Mass');



%---------- Altimeter range --------------------------------------
range_ice_sheet_20Hz_ku_name = 'range_ice_sheet_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,range_ice_sheet_20Hz_ku_name,int32_type,ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Corrected "ice sheet" altimeter range at 20 Hz for Ku band');
netcdf.putAtt(L2.ncid,id_aux,units_att,meters_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'LRM/PLRM modes: ice sheet (CFI) retracking, SAR mode: ice sheet margin retracking. Instrumental corrections included. Geophysical corrections included too.');

range_ice_sheet_20Hz_ku_mean_error_name = 'range_ice_sheet_20Hz_ku_mean_error';
id_aux = netcdf.defVar(L2.ncid,range_ice_sheet_20Hz_ku_mean_error_name,int32_type,single_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Mean error of the difference between the original values and a fitting of these');
netcdf.putAtt(L2.ncid,id_aux,units_att,meters_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_RANGE),'samples.'])

range_ice_sheet_20Hz_ku_rmse_name = 'range_ice_sheet_20Hz_ku_rmse';
id_aux = netcdf.defVar(L2.ncid,range_ice_sheet_20Hz_ku_rmse_name,int32_type,single_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Root Mean Squared Error of the difference between the original values and a fitting of these');
netcdf.putAtt(L2.ncid,id_aux,units_att,meters_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_RANGE),'samples.'])



%---------- Sigma0 --------------------------------
sig0_ice_sheet_20Hz_ku_name = 'sig0_ice_sheet_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,sig0_ice_sheet_20Hz_ku_name,int16_type,ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Backscatter coefficient sigma0');
netcdf.putAtt(L2.ncid,id_aux,units_att,dB_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'LRM/PLRM modes : ice sheet (CFI) retracking, SAR mode : ice sheet margin retracking. Instrumental corrections included : AGC instrumental errors correction (agc_cor_[x1]_[x2]) and internal calibration correction (sig0_cal_[x1]_[x2])');

sig0_ice_sheet_20Hz_ku_mean_error_name = 'sig0_ice_sheet_20Hz_ku_mean_error';
id_aux = netcdf.defVar(L2.ncid,sig0_ice_sheet_20Hz_ku_mean_error_name,int32_type,single_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Mean error of the difference between the original values and a fitting of these');
netcdf.putAtt(L2.ncid,id_aux,units_att,dB_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_SIG0),'samples.'])

sig0_ice_sheet_20Hz_ku_rmse_name = 'sig0_ice_sheet_20Hz_ku_rmse';
id_aux = netcdf.defVar(L2.ncid,sig0_ice_sheet_20Hz_ku_rmse_name,int32_type,single_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Root Mean Squared Error of the difference between the original values and a fitting of these');
netcdf.putAtt(L2.ncid,id_aux,units_att,dB_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_SIG0),'samples.'])



%----------  Elevation of echoing points --------------------------------
elevation_ice_sheet_20Hz_ku_name = 'elevation_ice_sheet_20Hz_ku';
id_aux = netcdf.defVar(L2.ncid,elevation_ice_sheet_20Hz_ku_name,int32_type,ku_rec_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Elevation of echoing points. Instrumental corrections included. Geophysical corrections included too.');
netcdf.putAtt(L2.ncid,id_aux,units_att,meters_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,'LRM/PLRM modes : ice sheet (CFI) retracking, SAR mode : ice sheet margin retracking. Instrumental corrections included : AGC instrumental errors correction (agc_cor_[x1]_[x2]) and internal calibration correction (sig0_cal_[x1]_[x2])');

elevation_ice_sheet_20Hz_ku_mean_error_name = 'elevation_ice_sheet_20Hz_ku_mean_error';
id_aux = netcdf.defVar(L2.ncid,elevation_ice_sheet_20Hz_ku_mean_error_name,int32_type,single_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Mean error of the difference between the original values and a fitting of these');
netcdf.putAtt(L2.ncid,id_aux,units_att,meters_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_ELEV),'samples.'])

elevation_ice_sheet_20Hz_ku_rmse_name = 'elevation_ice_sheet_20Hz_ku_rmse';
id_aux = netcdf.defVar(L2.ncid,elevation_ice_sheet_20Hz_ku_rmse_name,int32_type,single_dimension);
netcdf.putAtt(L2.ncid,id_aux,long_name_att,'Root Mean Squared Error of the difference between the original values and a fitting of these');
netcdf.putAtt(L2.ncid,id_aux,units_att,meters_units);
netcdf.putAtt(L2.ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(L2.ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(L2.ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_ELEV),'samples.'])



%----------  Global Attributes definition -----------------------------------
%---- attributes inherited from Sentinel-3 product description-------------
id_aux = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(L2.ncid,id_aux,'creation_time',date_creation);
netcdf.putAtt(L2.ncid,id_aux,'Conventions',netcdf_v4_format);
netcdf.putAtt(L2.ncid,id_aux,'altimeter_sensor_name',altimeter_sensor_name);
netcdf.putAtt(L2.ncid,id_aux,'first_meas_time',first_meas_time);
netcdf.putAtt(L2.ncid,id_aux,'last_meas_time',last_meas_time);
netcdf.putAtt(L2.ncid,id_aux,'first_meas_lat',first_meas_lat);
netcdf.putAtt(L2.ncid,id_aux,'last_meas_lat',last_meas_lat);
netcdf.putAtt(L2.ncid,id_aux,'first_meas_lon',first_meas_lon);
netcdf.putAtt(L2.ncid,id_aux,'last_meas_lon',last_meas_lon);
netcdf.putAtt(L2.ncid,id_aux,'semi_major_ellipsoid_axis',semi_major_ellipsoid_axis);
netcdf.putAtt(L2.ncid,id_aux,'ellipsoid_flattening',ellipsoid_flattening);
%--------------- add the attributes related to intermediate product--------
netcdf.putAtt(L2.ncid,id_aux,'orbit_cycle_num',orbit_cycle_num);
netcdf.putAtt(L2.ncid,id_aux,'orbit_REL_Orbit',orbit_REL_Orbit);
netcdf.putAtt(L2.ncid,id_aux,'orbit_ABS_Orbit_Start',orbit_ABS_Orbit_Start);

netcdf.endDef(L2.ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);

end