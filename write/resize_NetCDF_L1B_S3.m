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
% 2.1 Changed N_samples_sar_chd by N_samples
% 2.2 added cnf flag processing_mode_cnf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function resize_NetCDF_L1B_S3(files,N_total_surfaces)
global c_cst bw_ku_chd  zp_fact_range_cnf sec_in_day_cst pi_cst
global mission mode N_samples N_max_beams_stack_chd processing_mode_cnf
global compute_L1_stadistics_cnf include_wfms_aligned

global ACDC_application_cnf
ncid = netcdf.open(files.filename_netCDF,'WRITE');
ncid_2remove=netcdf.open(files.filename_netCDF_2remove,'NOWRITE');


%% PACKING L1B

%----------A. Time variables ----------------------------------------------
var_id = netcdf.inqVarID(ncid,'time_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'UTC_day_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'UTC_sec_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

%----------B. Orbit and attitude variables ------------------------------
var_id = netcdf.inqVarID(ncid,'lat_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'lon_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'alt_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'orb_alt_rate_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'satellite_mispointing_l1b_sar_echo_ku');
aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[3 N_total_surfaces]);
netcdf.putVar(ncid,var_id,[0 0],[3 N_total_surfaces],aux); clear aux;



%----------C. Flag time variables --------------------------------------

%----------D. Position/Velocity variables ------------------------------
var_id = netcdf.inqVarID(ncid,'x_pos_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'y_pos_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'z_pos_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'x_vel_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'y_vel_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'z_vel_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

%----------E. Navigation Bulletin ------------------------------

% var_id = netcdf.inqVarID(ncid,'seq_count_l1b_echo_sar_ku');
% netcdf.putVar(ncid,var_id,i_surf,seq_count_l1b_echo_sar_ku);

%----------F. Operating instrument & tracking --------------------------
var_id = netcdf.inqVarID(ncid,'oper_instr_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'SAR_mode_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

%---------- G. H0, COR2 and AGC ----------------------------------------
var_id = netcdf.inqVarID(ncid,'h0_applied_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
% aux=ncread(files.filename_netCDF_2remove,'h0_applied_l1b_echo_sar_ku',1,N_total_surfaces);
% ncwrite(files.filename_netCDF,'h0_applied_l1b_echo_sar_ku',aux,1); clear aux;

var_id = netcdf.inqVarID(ncid,'cor2_applied_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'agccode_ku_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

%---------- H. Surface type -----------------------------------------------
var_id = netcdf.inqVarID(ncid,'surf_type_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

%---------- I. Altimeter range and Corrections ----------------------------
var_id = netcdf.inqVarID(ncid,'range_ku_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

if include_wfms_aligned
    var_id = netcdf.inqVarID(ncid,'range_ku_surf_aligned_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

var_id = netcdf.inqVarID(ncid,'uso_cor_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'int_path_cor_ku_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(mode,'SIN') && strcmp(mission,'CR2') && strcmp(processing_mode_cnf,'SIN')
    var_id = netcdf.inqVarID(ncid,'int_path_2_cor_ku_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
end

var_id = netcdf.inqVarID(ncid,'range_rate_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

%---------- J. AGC and Sigma0 scalings --------------------------------
var_id = netcdf.inqVarID(ncid,'scale_factor_ku_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

%---------- K. Stack characterization--------------------------------------
var_id = netcdf.inqVarID(ncid,'nb_stack_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'nb_stack_start_stop_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'look_angle_start_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'look_angle_stop_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'doppler_angle_start_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
 
var_id = netcdf.inqVarID(ncid,'doppler_angle_stop_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'pointing_angle_start_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
 
var_id = netcdf.inqVarID(ncid,'pointing_angle_stop_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

if compute_L1_stadistics_cnf
    var_id = netcdf.inqVarID(ncid,'skew_stack_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'kurt_stack_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'stdev_stack_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
%     aux=ncread(files.filename_netCDF_2remove,'stdev_stack_l1b_echo_sar_ku',1,N_total_surfaces);
%     ncwrite(files.filename_netCDF,'stdev_stack_l1b_echo_sar_ku',aux,1); clear aux;
    
    
    var_id = netcdf.inqVarID(ncid,'gaussian_fitting_centre_look_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'gaussian_fitting_centre_pointing_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'beam_form_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
%     aux=ncread(files.filename_netCDF_2remove,'beam_form_l1b_echo_sar_ku',1,N_total_surfaces);
%     ncwrite(files.filename_netCDF,'beam_form_l1b_echo_sar_ku',aux,1); clear aux;
    
end

%------------- L. Altimeter engineering variables -------------------------
var_id = netcdf.inqVarID(ncid,'altimeter_clock_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'pri_lrm_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
% aux=ncread(files.filename_netCDF_2remove,'pri_lrm_l1b_echo_sar_ku',1,N_total_surfaces);
% ncwrite(files.filename_netCDF,'pri_lrm_l1b_echo_sar_ku',aux,1); clear aux;


%------------- M. Waveform related variables -----------------------------
var_id = netcdf.inqVarID(ncid,'i2q2_meas_ku_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
% aux=ncread(files.filename_netCDF_2remove,'i2q2_meas_ku_l1b_echo_sar_ku',[1 1],[N_samples*zp_fact_range_cnf N_total_surfaces]);
% ncwrite(files.filename_netCDF,'i2q2_meas_ku_l1b_echo_sar_ku',aux,[1 1]); clear aux;

if include_wfms_aligned
    var_id = netcdf.inqVarID(ncid,'i2q2_meas_ku_wdcorr_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
    netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
%     aux=ncread(files.filename_netCDF_2remove,'i2q2_meas_ku_wdcorr_l1b_echo_sar_ku',[1 1],[N_samples*zp_fact_range_cnf N_total_surfaces]);
%     ncwrite(files.filename_netCDF,'i2q2_meas_ku_wdcorr_l1b_echo_sar_ku',aux,[1 1]); clear aux;
end

var_id = netcdf.inqVarID(ncid,'stack_mask_range_bin_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_max_beams_stack_chd N_total_surfaces]);
netcdf.putVar(ncid,var_id,[0 0],[N_max_beams_stack_chd N_total_surfaces],aux); clear aux;


var_id = netcdf.inqVarID(ncid,'waveform_scale_factor_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;


if strcmp(mode,'SIN') && strcmp(mission,'CR2') && strcmp(processing_mode_cnf,'SIN') 
    var_id = netcdf.inqVarID(ncid,'i2q2_meas_ku_l1b_echo_sar_ku_2');
    aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
    netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'coherence_meas_ku_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
    netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'phase_difference_meas_ku_l1b_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
    netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'waveform_scale_factor_l1b_echo_sar_ku_2');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
             
    
    if include_wfms_aligned
        var_id = netcdf.inqVarID(ncid,'i2q2_meas_ku_wdcorr_l1b_echo_sar_ku_2');
        aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
        netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
        
        var_id = netcdf.inqVarID(ncid,'coherence_meas_ku_wdcorr_l1b_echo_sar_ku');
        aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
        netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
        
        var_id = netcdf.inqVarID(ncid,'phase_difference_meas_ku_wdcorr_l1b_echo_sar_ku');
        aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
        netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
    end
           
end

%added EM 04.10.2016
% %------------- P. ACDC variables ----------------------------------------
if ACDC_application_cnf
    %-------------- Waveform -----------------------------------------------    
    %ACDC waveform
    var_id = netcdf.inqVarID(ncid,'i2q2_meas_ku_ACDC_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
    netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
    
    %fitted ACDC waveform    
    var_id = netcdf.inqVarID(ncid,'i2q2_fit_ku_ACDC_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
    netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
    
    
    %fitted L1B waveform
    var_id = netcdf.inqVarID(ncid,'i2q2_fit_ku_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
    netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
    
    
    %-------------- Scaling -----------------------------------------------
    %ACDC waveform
    var_id = netcdf.inqVarID(ncid,'waveform_scale_factor_ACDC_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    %fitted ACDC waveform
    var_id = netcdf.inqVarID(ncid,'waveform_scale_factor_fit_ACDC_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    %fitted L1B waveform
    var_id = netcdf.inqVarID(ncid,'waveform_scale_factor_fit_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    
    
    %-------------- Range index -------------------------------------------
    var_id = netcdf.inqVarID(ncid,'range_index_ACDC_echo_sar_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces]);
    netcdf.putVar(ncid,var_id,[0 0],[N_samples*zp_fact_range_cnf N_total_surfaces],aux); clear aux;
    
    %-------------- Geophysical retrievals --------------------------------
    %******************* retracked range **********************************
    var_id = netcdf.inqVarID(ncid,'retracked_range_ACDC_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    %******************* epoch ********************************************
    var_id = netcdf.inqVarID(ncid,'epoch_ACDC_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    %******************** SWH  ********************************************
    var_id = netcdf.inqVarID(ncid,'swh_ACDC_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    var_id = netcdf.inqVarID(ncid,'swh_ACDC_PRE_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    %******************** SSH  ********************************************
    var_id = netcdf.inqVarID(ncid,'ssh_ACDC_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    var_id = netcdf.inqVarID(ncid,'ssh_ACDC_PRE_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    %******************** SIGMA0 ******************************************    
    var_id = netcdf.inqVarID(ncid,'sig0_ACDC_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'sig0_ACDC_PRE_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    %********************* Pu: fitted peak power **************************    
    var_id = netcdf.inqVarID(ncid,'Pu_analytical_ACDC_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    %********************** Pearson Correlation Coefficient ***************
    var_id = netcdf.inqVarID(ncid,'Pearson_corr_ACDC_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(ncid,'Pearson_corr_ACDC_PRE_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;    
    %************************ Fitting flag ********************************    
    var_id = netcdf.inqVarID(ncid,'Flag_fitting_ACDC_20_ku');
    aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;
      
end


%------------- N. Geophysical Corrections variables -----------------------

var_id = netcdf.inqVarID(ncid,'dry_tropo_correction_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'wet_tropo_correction_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'inverse_baro_correction_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'Dynamic_atmospheric_correction_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'GIM_iono_correction_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'model_iono_correction_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'ocean_equilibrium_tide_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'long_period_tide_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'ocean_loading_tide_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'solid_earth_tide_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'geocentric_polar_tide_l1b_echo_sar_ku');
aux=netcdf.getVar(ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(ncid,var_id,0,N_total_surfaces,aux); clear aux;

netcdf.close(ncid_2remove);
netcdf.close(ncid);