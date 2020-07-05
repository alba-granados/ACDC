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
% 2.1 added cnf flag processing_mode_cnf
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function resize_NetCDF_L1A_S3(filename,start, count, files) % alba

global N_samples 

% ncid = netcdf.open('./results/data/measurement_l1a_reduced.nc','WRITE');
ncid = netcdf.open(strcat(strcat(files.resultPath,'data/'),'measurement_l1a_reduced.nc'),'WRITE'); % alba
ncid_L1A_big = netcdf.open(filename,'NOWRITE');


var_id_reduced = netcdf.inqVarID(ncid,'seq_count_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'seq_count_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,int16(aux)); clear aux;


var_id_reduced = netcdf.inqVarID(ncid,'UTC_day_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'UTC_day_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'UTC_sec_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'UTC_sec_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'uso_cor_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'uso_cor_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'burst_count_prod_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'burst_count_prod_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'burst_count_prod_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'burst_count_prod_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'pitch_sat_pointing_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'pitch_sat_pointing_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'yaw_sat_pointing_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'seq_count_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,int16(aux)); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'cog_cor_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'cog_cor_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'x_pos_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'x_pos_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'y_pos_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'y_pos_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'z_pos_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'z_pos_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'x_vel_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'x_vel_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'y_vel_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'y_vel_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'z_vel_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'z_vel_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'alt_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'alt_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'orb_alt_rate_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'orb_alt_rate_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'roll_sral_mispointing_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'roll_sral_mispointing_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'pitch_sral_mispointing_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'pitch_sral_mispointing_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'yaw_sral_mispointing_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'yaw_sral_mispointing_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'cor2_applied_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'cor2_applied_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'h0_applied_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'h0_applied_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,int32(aux)); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'burst_count_cycle_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'burst_count_cycle_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'int_path_cor_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'int_path_cor_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'uso_cor_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'uso_cor_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;

var_id_reduced = netcdf.inqVarID(ncid,'range_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'range_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'h0_nav_dem_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'h0_nav_dem_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,int32(aux)); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'cor2_nav_dem_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'cor2_nav_dem_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'agc_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'agc_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'sig0_cal_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'sig0_cal_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'scale_factor_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'scale_factor_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;


var_id_reduced = netcdf.inqVarID(ncid,'burst_power_cor_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'burst_power_cor_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,[0 start], [64 count]);
netcdf.putVar(ncid,var_id_reduced,[0 0],[64 count],int32(aux)); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'burst_phase_cor_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'burst_phase_cor_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,[0 start], [64 count]);
netcdf.putVar(ncid,var_id_reduced,[0 0],[64 count],aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'gprw_meas_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'gprw_meas_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,[0 0 start], [N_samples 3 count]);
netcdf.putVar(ncid,var_id_reduced,[0 0 0],[N_samples 3 count],int32(aux)); clear aux;


var_id_reduced = netcdf.inqVarID(ncid,'i_meas_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'i_meas_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,[0 0 start], [N_samples 64 count]);
netcdf.putVar(ncid,var_id_reduced,[0 0 0],[N_samples 64 count],aux); clear aux;
var_id_reduced = netcdf.inqVarID(ncid,'q_meas_ku_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'q_meas_ku_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,[0 0 start], [N_samples 64 count]);
netcdf.putVar(ncid,var_id_reduced,[0 0 0],[N_samples 64 count],aux); clear aux;


var_id_reduced = netcdf.inqVarID(ncid,'time_l1a_echo_sar_ku');
var_id_big = netcdf.inqVarID(ncid_L1A_big,'time_l1a_echo_sar_ku');
aux=netcdf.getVar(ncid_L1A_big,var_id_big,start,count);
netcdf.putVar(ncid,var_id_reduced,0,count,aux); clear aux;


netcdf.close(ncid_L1A_big);
netcdf.close(ncid);