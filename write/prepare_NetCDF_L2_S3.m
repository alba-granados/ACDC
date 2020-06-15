%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-1B product using Sentinel-3 like format
% Ref: Product Data Format Specification - SRAL/MWR Level 1 & 2 Instrument
% Products issue 2.0
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
% 2.1 changed N_samples_sar_chd by N_samples
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function prepare_NetCDF_L2_S3(L2)

%global netcdf_type
% ncid = netcdf.open(L2.filename_netCDF,'WRITE');

% ku_rec_dimension = netcdf.defDim(ncid,'time_l1b_echo_sar_ku',N_surfs_loc_estimated);
% nl_dimension = netcdf.defDim(ncid,'Nl',N_max_beams_stack_chd);
% ns_dimension = netcdf.defDim(ncid,'Ns',N_samples*zp_fact_range_cnf);
% space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);

% t6 = tic;

% date_creation = datestr(now, '_yyyymmdd_');
% switch mission
%     case 'CR2'        
%         switch mode
%             case 'SAR'                
%                 files.filename_netCDF = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',...
%                                 files.sph.product_info.product_id(20:20+30),...
%                                 date_creation,...
%                                 'isd','.nc');
%         end
%     case {'S3_','S3A','S3B'} 
%     case {'S6_'} 
% end


% dimensions_key = 'Dimensions';
% format_key = 'Format';
% data_type_key = 'DataType';
% fill_value_key='FillValue';
% 
% netcdf_v4_format = 'netcdf4';
% 
% ku_rec_dimension = 'time_l1b_echo_sar_ku');
% nl_dimension = 'Nl';
% nl_dimension_size = N_max_beams_stack_chd;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;



% int8_type = 'int8';
% uint8_type = 'uint8';
% int16_type = 'int16';
% uint16_type = 'uint16';
% int32_type = 'int32';
% uint32_type = 'uint32';
% uint64_type = 'uint64';
% float_type = 'single';
% double_type='double';


%% CODING L2

%----------A. Time variables ----------------------------------------------
time_20Hz_ku    = L2.ISD_time_surf(L2.idx_int_ISD); %leap seconds are already taken into account in L1B time values
UTC_day_20Hz_ku = L2.ISD_UTC_day(L2.idx_int_ISD);
UTC_sec_20Hz_ku = L2.ISD_UTC_sec(L2.idx_int_ISD);

% UTC_day_20Hz_ku = int16(floor(L1BS.time_surf./ sec_in_day_cst));
% UTC_sec_20Hz_ku = double((time_20Hz_ku - double(UTC_day_20Hz_ku) * sec_in_day_cst));


%---------- Orbit variables ------------------------------
lat_20Hz_ku = int32(L2.ISD_lat_surf(L2.idx_int_ISD)         * 1e6);
lon_20Hz_ku = int32(L2.ISD_lon_surf(L2.idx_int_ISD)         * 1e6); 
alt_20Hz_ku = int32((L2.ISD_alt_sat(L2.idx_int_ISD)-700000) * 1e4);


%---------- Altimeter range -------------------------
range_ice_sheet_20Hz_ku                 = int32((L2.ISD_range-700000)   * 1e4);
range_ice_sheet_20Hz_ku_mean_error      = int32(L2.range_mean_error     * 1e4);
range_ice_sheet_20Hz_ku_rmse            = int32(L2.range_rmse           * 1e4);

%---------- Sigma0 --------------------------------
sig0_ice_sheet_20Hz_ku                  = int16(L2.ISD_sigma0           * 1e2);
sig0_ice_sheet_20Hz_ku_mean_error       = int16(L2.sig0_mean_error      * 1e2);
sig0_ice_sheet_20Hz_ku_rmse             = int16(L2.sig0_rmse            * 1e2);

%---------- Elevation of echoing points -------------------------
elevation_ice_sheet_20Hz_ku             = int32((L2.ISD_SSH)            * 1e4);
elevation_ice_sheet_20Hz_ku_mean_error  = int32(L2.elev_mean_error      * 1e4);
elevation_ice_sheet_20Hz_ku_rmse        = int32(L2.elev_rmse            * 1e4);



%% PACKING L2
%----------A. Time variables ----------------------------------------------
var_id = netcdf.inqVarID(L2.ncid,'time_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,time_20Hz_ku);

var_id = netcdf.inqVarID(L2.ncid,'UTC_day_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,UTC_day_20Hz_ku);

var_id = netcdf.inqVarID(L2.ncid,'UTC_sec_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,UTC_sec_20Hz_ku);


%------------------------ Orbit variables ---------------------------------
var_id = netcdf.inqVarID(L2.ncid,'lat_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,lat_20Hz_ku);

var_id = netcdf.inqVarID(L2.ncid,'lon_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,lon_20Hz_ku);

var_id = netcdf.inqVarID(L2.ncid,'alt_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,alt_20Hz_ku);


%------------------------ Altimeter range ---------------------------------
var_id = netcdf.inqVarID(L2.ncid,'range_ice_sheet_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,range_ice_sheet_20Hz_ku);

var_id = netcdf.inqVarID(L2.ncid,'range_ice_sheet_20Hz_ku_mean_error');
netcdf.putVar(L2.ncid,var_id,range_ice_sheet_20Hz_ku_mean_error);

var_id = netcdf.inqVarID(L2.ncid,'range_ice_sheet_20Hz_ku_rmse');
netcdf.putVar(L2.ncid,var_id,range_ice_sheet_20Hz_ku_rmse);


%--------------------- Sigma0 scaling factor ------------------------------
var_id = netcdf.inqVarID(L2.ncid,'sig0_ice_sheet_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,sig0_ice_sheet_20Hz_ku);

var_id = netcdf.inqVarID(L2.ncid,'sig0_ice_sheet_20Hz_ku_mean_error');
netcdf.putVar(L2.ncid,var_id,sig0_ice_sheet_20Hz_ku_mean_error);

var_id = netcdf.inqVarID(L2.ncid,'sig0_ice_sheet_20Hz_ku_rmse');
netcdf.putVar(L2.ncid,var_id,sig0_ice_sheet_20Hz_ku_rmse);


%----------------  Elevation of echoing points ----------------------------
var_id = netcdf.inqVarID(L2.ncid,'elevation_ice_sheet_20Hz_ku');
netcdf.putVar(L2.ncid,var_id,elevation_ice_sheet_20Hz_ku);

var_id = netcdf.inqVarID(L2.ncid,'elevation_ice_sheet_20Hz_ku_mean_error');
netcdf.putVar(L2.ncid,var_id,elevation_ice_sheet_20Hz_ku_mean_error);

var_id = netcdf.inqVarID(L2.ncid,'elevation_ice_sheet_20Hz_ku_rmse');
netcdf.putVar(L2.ncid,var_id,elevation_ice_sheet_20Hz_ku_rmse);



% %----------  Global Attributes definition -----------------------------------
% %---- attributes inherited from Sentinel-3 product description-------------
% 
% ncwriteatt(L2.filename_netCDF,'/','Conventions',netcdf_v4_format);
% ncwriteatt(L2.filename_netCDF,'/','altimeter_sensor_name',altimeter_sensor_name);
% ncwriteatt(L2.filename_netCDF,'/','first_meas_time',first_meas_time);
% ncwriteatt(L2.filename_netCDF,'/','last_meas_time',last_meas_time);
% ncwriteatt(L2.filename_netCDF,'/','first_meas_lat',first_meas_lat);
% ncwriteatt(L2.filename_netCDF,'/','last_meas_lat',last_meas_lat);
% ncwriteatt(L2.filename_netCDF,'/','first_meas_lon',first_meas_lon);
% ncwriteatt(L2.filename_netCDF,'/','last_meas_lon',last_meas_lon);
% ncwriteatt(L2.filename_netCDF,'/','semi_major_ellipsoid_axis',semi_major_ellipsoid_axis);
% ncwriteatt(L2.filename_netCDF,'/','ellipsoid_flattening',ellipsoid_flattening);
% %--------------- add the attributes related to intermediate product--------
% ncwriteatt(L2.filename_netCDF,'/','orbit_cycle_num',orbit_cycle_num);
% ncwriteatt(L2.filename_netCDF,'/','orbit_REL_Orbit',orbit_REL_Orbit);
% ncwriteatt(L2.filename_netCDF,'/','orbit_ABS_Orbit_Start',orbit_ABS_Orbit_Start);


netcdf.close(L2.ncid);

end