% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% CR2 SARIn 
% ---------------------------------------------------------
% Objective: Define Configuration parameters
% 
% INPUTs : - 
% OUTPUTs: -
%
% ----------------------------------------------------------
% Author:    Roger Escolà  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Mònica Roca   / isardSAT
% Last rev.: Mònica Roca   / isardSAT (11/05/2015)
%
%
% Version  record
% 1.0 2015/01/01 Imported code from S6
% 1.1 2016/02/18 Added hamming_window_cnf and zp_fact_azimut_cnf
% 1.2 2016/03/24 Added mask_flag
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filtering 
global mask_flag
mask_flag = 0;

%% Processing and testing
global clock_unit_conv_cnf cor2_bit_shift_cnf  
global sigma_alt_surf_th_cnf


%% PRELIMINARY DATATION

%% PRELIMINARY WINDOW DELAY

%% FINAL BURST DATATION

%% ONBOARD REVERSION
global height_rate_application_cnf FAI_application_cnf
height_rate_application_cnf=0;
FAI_application_cnf=1;

%% INSTRUMENT GAIN

%% WAVEFORMS CORRECTION
global CAL2_flag_cnf CAL1p2p_flag_cnf
CAL2_flag_cnf=1;
CAL1p2p_flag_cnf=1;

%% SURFACE LOCATION

global N_surf_samp_interpol_cnf 
global N_points_spline_surf_cnf N_points_spline_orbit_cnf
global smooth_fact_surf_pos_cnf 
global smooth_fact_orbit_pos_cnf smooth_fact_orbit_vel_cnf
global accu_orbit_alt_cnf accu_orbit_lat_cnf
global sigma_alt_surf_geoloc_th_cnf

N_surf_samp_interpol_cnf = 20;
N_points_spline_surf_cnf = 10;
smooth_fact_surf_pos_cnf = 0;
N_points_spline_orbit_cnf = N_points_spline_surf_cnf;
smooth_fact_orbit_pos_cnf = 0;
smooth_fact_orbit_vel_cnf = 0;
accu_orbit_alt_cnf = 0;
accu_orbit_lat_cnf = 0.0001;
sigma_alt_surf_geoloc_th_cnf = 100 * N_surf_samp_interpol_cnf; %SURFACE LOCATIONS (HIGH OR LOW ROUGHNESS)

global create_surfaces_backward_cnf create_surfaces_foreward_cnf
global backward_margin_cnf foreward_margin_cnf
create_surfaces_backward_cnf = 0;
create_surfaces_foreward_cnf = 0;
backward_margin_cnf = 32+16;
foreward_margin_cnf = 32+16;

global trp_flag_cnf trp_flag_method_cnf
global lat_trp lon_trp alt_trp

trp_flag_cnf            = 0; % TRP PROCESSING
trp_flag_method_cnf     = 1;

if (trp_flag_cnf)
    % Crete TRP
     lat_trp = 35.3379302808;
     lon_trp = 23.7795182869;
     alt_trp = 1048.8184;
    % Svalbard TRP
%      lat_trp = 78.23052306;
%      lon_trp = 15.39376997;
%      alt_trp = 492.772;
end

%% BEAM ANGLES

%% AZIMUTH PROCESSING
global hamming_window_cnf gain_scale_hamm
global zp_fact_azimut_cnf zp_fact_range_cnf 
global force_exact_method_cnf hamm_win
global N_ku_pulses_burst_chd

sigma_alt_surf_th_cnf = 1000; %AZIMUTH PROCESSING (HIGH OR LOW VARIABILITY)
zp_fact_azimut_cnf = 1;
hamming_window_cnf = 1;
force_exact_method_cnf=0;
if(hamming_window_cnf)
    %compute the gain scaling of the hamming window as defined by dinardo
%     hamm = hamming(N_ku_pulses_burst_chd+1); %we want the hamming centered in sample #0 (if the samples go from -32 to 31)
%     hamm_win = hamm(1:N_ku_pulses_burst_chd);
%     clear hamm;
    %Hamming according to Dinardo
    c1=0.08;
    c2=0.92;
    hamm_win=(c1+c2.*(cos(pi.*(0:1:N_ku_pulses_burst_chd-1)./(N_ku_pulses_burst_chd-1)-pi/2)).^2).';
%     %Unitary energy window
    %gain_scale_hamm=-10*log10(sum(hamm_win.^2));
    %Dinardo option
    gain_scale_hamm=20*log10(mean(hamm_win));
    
else
    gain_scale_hamm=0;
    hamm_win=[];
end

%% GEOMETRY CORRECTION

global win_delay_ref_cnf elevation_ref_cnf
global apply_doppler apply_slant apply_wd
win_delay_ref_cnf = 0;       %0 = normal alignment
                             %1 = optimal alignment for ICE (wrt max beam accumulated power) CR2
                             %2 = force the alignment to the closest of a given window_delay: Lakes
                             %3 = force to the first one.
                             %4 = force to the minimum one.
                             if(win_delay_ref_cnf==2)
                                 elevation_ref_cnf = 1048.8184; %force the window delay to a reference elevation
                             end
apply_doppler = 1;
apply_slant = 1;
apply_wd = 1;

%% RANGE COMPRESSION
global coherence_threshold_cnf 
global window_rg_cnf window_range_SAR window_range_SARin gain_scale_win_rg_SAR gain_scale_win_rg_SARin
global N_samples_sar_chd N_samples_sin_chd 
global stack_interferometry_flag_cnf
coherence_threshold_cnf = 10;
stack_interferometry_flag_cnf = 1;
zp_fact_range_cnf = 1;
window_rg_cnf = 0;

global extend_window_cnf
extend_window_cnf = 0;
if(window_rg_cnf)
    %compute the gain scaling of the hamming window as defined by dinardo
%     hamm = hamming(N_ku_pulses_burst_chd+1); %we want the hamming centered in sample #0 (if the samples go from -32 to 31)
%     hamm_win = hamm(1:N_ku_pulses_burst_chd);
%     clear hamm;
    %Hamming according to Dinardo
    c1=0.08;
    c2=0.92;
    window_range_SAR=(c1+c2.*(cos(pi.*(0:1:N_samples_sar_chd-1)./(N_samples_sar_chd-1)-pi/2)).^2);
    window_range_SARin=(c1+c2.*(cos(pi.*(0:1:N_samples_sin_chd-1)./(N_samples_sin_chd-1)-pi/2)).^2);
%     %Unitary energy window
    %gain_scale_hamm=-10*log10(sum(hamm_win.^2));
    %Dinardo option
    gain_scale_win_rg_SAR=20*log10(mean(window_range_SAR));
    gain_scale_win_rg_SARin=20*log10(mean(window_range_SARin));
    
else
    gain_scale_win_rg_SAR=0;
    gain_scale_win_rg_SARin=0;
    window_range_SAR=[];
    window_range_SARin=[];
end

%% MULTILOOKING
global compute_L1_stadistics_cnf N_beams_sub_stack_cnf
global range_index1_cnf range_index2_cnf
global k_noise_top_cnf k_noise_floor_cnf
global use_zeros_cnf delete_ambiguities_cnf
global coherency_mask_cnf  avoid_wrapps_cnf apply_stack_mask_cnf
global avoid_beams_mask_allzeros_cnf
global mask_look_angles_CR2_cnf mask_look_angles_CR2_min_chd mask_look_angles_CR2_max_chd
global avoid_noisy_beams_cnf method_noise_est_cnf noise_estimation_window param_std_noise_thres_cnf

compute_L1_stadistics_cnf = 0;
N_beams_sub_stack_cnf = 5;
range_index1_cnf = 12;
range_index2_cnf = 16;
k_noise_top_cnf = 7;
k_noise_floor_cnf = 0;
use_zeros_cnf = 1;
apply_stack_mask_cnf=1;
delete_ambiguities_cnf=0;
coherency_mask_cnf = 1;

avoid_wrapps_cnf = 0;

avoid_beams_mask_allzeros_cnf=1;

mask_look_angles_CR2_max_chd=0.6*pi/180;
mask_look_angles_CR2_min_chd=-0.6*pi/180;
mask_look_angles_CR2_cnf=0;

avoid_noisy_beams_cnf=0; % waveforms of outer noisy beams in the stack are discarded 
                         % following Sentinel-3 baseline: threshold=noise_mean+3*std_noise       
param_std_noise_thres_cnf=3.0; %number of times the standard deviation of the noise used in the thresholding
noise_estimation_window=[12 16]; %first and last non-zero padded noise samples
method_noise_est_cnf='after_geom_corr'; % String: 
                       % 'after_geom_corr': noise estimation is performed
                       % right after the geometry corrections as done in
                       % Sentinel-3
                       % 'before_geom_corr': noise estimation is performed
                       % right before the geometry corrections as done in
                       % Sentinel-3                      
% method_noise_est_cnf='Beam_based'; % String: 
%                        % 'Beam_based': for each beam its related noise
%                        % statistics are estimated and accordingly applied
%                        % 'Stack_based': noise statistics estimated using
%                        % all the samples over the stack for the defined
%                        % noise estimation window --> thresholding for
%                        % different beams based on a single mean and std
% method_noisy_thresholding_cnf='Average_based'; % String:
%                                 % 'Peak_based': Compare peak with threshold
%                                 % 'Average_based': Compare average with
%                                 % threshold





            

%% not used?
clock_unit_conv_cnf = 256;
cor2_bit_shift_cnf = 4;

%% Product writing
global include_wfms_aligned
include_wfms_aligned=0;
% global netcdf_type
% netcdf_type='netcdf4'; %string indicating the type of netcdf 3 or 4: netcdf3 or netcdf4: 
% %Consequence on the uint not available in the netcdf3 definition need to
% %enlarge to int ouf double of bits to accomodate it --> enlarge the product

%% Verbosity: figures pop up or not
global figures_visible_cnf
figures_visible_cnf=0;

%% ------------ ACDC related parameters ------------------------------------
global ACDC_application_cnf 
ACDC_application_cnf=0;

if ACDC_application_cnf
    global cnf_p_ACDC
    %--------------------------------------------------------------------------
    %---------------------- PRELIMINARY ESTIMATION OF SWH ---------------------
    %--------------------------------------------------------------------------
    cnf_p_ACDC.method_pre_fitting=0; %Indicate which approach to be used for preliminary estimation of SWH
    %0: using whole stack (start to stop beams)
    %1: using region around specular beam
    %2: using region around maximum beam
    cnf_p_ACDC.num_left_beams_est_cnf= 10;%number of beams to the left for estimation when using methods 1 or 2;
    cnf_p_ACDC.num_right_beams_est_cnf= 10;%number of beams to the right for estimation when using methods 1 or 2;
    cnf_p_ACDC.preliminary_est_surf_downsampling=1; %every how much we need to re-estimate as a preliminary step the Hs using one of the above methods
    
    %--------------------------------------------------------------------------
    %----------------------- RETRACKING ML ACDC -------------------------------
    %--------------------------------------------------------------------------
    cnf_p_ACDC.weighting_win_type='gaussian'; % String to indicate the type of weighting 
                                    % window in the multilooking of the
                                    % re-order 1D AC stack: smoothing like
                                    % function
    cnf_p_ACDC.weighting_win_width=0.5; %width of the window in range bins
    
    
    %--------------------------------------------------------------------------
    %----------------------- L1B PROCESSING -----------------------------------
    %--------------------------------------------------------------------------
    cnf_p_ACDC.use_zeros_cnf=use_zeros_cnf;
    cnf_p_ACDC.ZP=zp_fact_range_cnf;  %to be modified to consider azimuth and range right now only range
    
    %--------------------------------------------------------------------------
    % ----------------------NOISE ESTIMATION ----------------------------------
    %--------------------------------------------------------------------------
    cnf_p_ACDC.Thn_flag          =   1;              % 1 - To account for ThN in estimation; 0 - Otherwise
    cnf_p_ACDC.Thn_w_first       =   1;%50;             % Gate number after IF Mask to start Thermal noise windowing; this is a subscript indice thus must be > 0
    cnf_p_ACDC.Thn_w_width       =   20;%20;             % Thermal noise window width in range bins
    cnf_p_ACDC.fit_noise         =   1; % 0 using the values in Thn_w_first and Thn_w_width or 1
    %to perform an adaptive estimation of the window width based
    %on pre-processing step estimation using the
    %multilooked waveform (derivative and leading edge detection)
    %PERFORMED WHENEVER PRE-PROCESSING FLAG IS
    %ACTIVE
    cnf_p_ACDC.threshold_noise   = 1e-3; %1e-3 %threshold used to estimate the samples used in the noise window estimation
    cnf_p_ACDC.max_iter_noise    = 100;
    cnf_p_ACDC.Thn_method        =  'ML'; % String:
    
    
    %--------------------------------------------------------------------------
    % ------------------- MODEL CONFIGURATION ---------------------------------
    %--------------------------------------------------------------------------
    % ------------------ Surface parameters -----------------------------------
    cnf_p_ACDC.rou_flag          =   0;              % NO TOQUIS If Gaussian PDF to be considered, then no need to fit roughness which for Ocean is equivalent to Skewness
    %if ~cnf_p_ACDC.rou_flag
    cnf_p_ACDC.rou           =   1e-2;           % Gaussian PDF show a roughness value in the order of 10-2
    cnf_p_ACDC.Hs            = 1e-6; % value of the Hs when fitting the MSS or roughness
    %end
    
    %------------------- Power waveform model ---------------------------------
    cnf_p_ACDC.power_wfm_model='complete'; %Define the model approximation of power wfm whether to compute
    % 'simple': Pkl=Bkl*sqrt(gl)*func_f0
    % 'complete':
    % Pkl=Bkl*sqrt(gl)*(func_f0+Tkl*gl*simga_s^2*func_f1);
    %-------------------- Multilooking option ---------------------------------
    cnf_p_ACDC.multilook_option='Cris_old'; %String:
    % -Cris_old: thermal noise is added to
    % each SL
    % -NormML_plus_noise: ML is normalized
    % and the thermal noise is added
    % outside
    % -------------------------------------------------------------------------
    % ---------------------- LOOK UP TABLES CONFIGURATION ---------------------
    % -------------------------------------------------------------------------
    % LUT tables in .mat should be in the processor folder
    cnf_p_ACDC.lut_flag          =   1;              % 0 or 1;if null do not use Look up tables
    cnf_p_ACDC.N_int             =   200; %not used            % number of intervals in the integration
    cnf_p_ACDC.lim_int           =   5.1; %not used           % the range of the integration
    cnf_p_ACDC.LUT_ximax         =   50; % maximum value of the gl*k parameter LUT
    cnf_p_ACDC.LUT_ximin         =   -10; % minimum value of the gl*k parameter LUT
    cnf_p_ACDC.LUT_step          =   1e-5; % step of gl*k parameter in the LUT
    
    %--------------------------------------------------------------------------
    %----------------- PRE-PROCESSING STAGE------------------------------------
    %--------------------------------------------------------------------------
    %provide better initial estimate of epoch--
    %include a pre_processing stage to compute the initial epoch if desired
    %based on a mid-leading edge estimation as a percentatge of the peak power
    %in the actual waveform and eventually include the estimation of the noise
    %floor adjusting automatically
    cnf_p_ACDC.pre_processing=1;
    cnf_p_ACDC.percent_leading_edge=87.0; % percentage of peak detect to establish the mid-point leading edge
    
    %----------------------------------------------------------------------
    %----------------------- ITERATIVE PROCESSING -------------------------
    %----------------------------------------------------------------------
    cnf_p_ACDC.num_iterations=1;
    
    % -------------------------------------------------------------------------
    % ------------- FITTING CONFIGURATION -------------------------------------
    %---- Initial values-------------------------------------------------------
    cnf_p_ACDC.ini_Epoch         =   35;             % 45 for BB, 145 for MEAS18
    cnf_p_ACDC.ini_Hs            =   2.0;
    cnf_p_ACDC.ini_Pu            =   1;
    cnf_p_ACDC.ini_rou          =   1e-2;
    cnf_p_ACDC.ini_error_epoch = 0.0; %differential error of the ACDC retracker
    cnf_p_ACDC.initial_param_fit_feedback_flag = 1; %activate the updating of the initial parameters per fitted waveform either using the sliding window or with acumulation of previous info
    cnf_p_ACDC.ini_Hs_rou_sliding_win_opt = 1;  % activate the averaging of Hs/roughness for next initial seed
    cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold = 3; %number of times std for removal (it is also used for comparison with the accumulated one)
    cnf_p_ACDC.ini_Hs_rou_sliding_win_size = 20;
    
    % -------------------SMOOTHING STACK AFTER AC -------------------------
    cnf_p_ACDC.stack_smooth_active=0;
    cnf_p_ACDC.stack_smooth_win=20;
    
    %------------------ Range indexation ----------------------------------
    cnf_p_ACDC.range_index_method='resampling'; % 'conventional': when creating the range vector used in fitting use directly 1:1:N_samples (even with ZP)
    % 'resampling': creating the
    % range vector as a resampling
    % of the original one (if there is a resampling)
    
    cnf_p_ACDC.fitting_fun_type='lsq'; % Define the fitting procedure used
    
    %String:
    %-'lsq': Using the lsqcurvfit based on levegender-Marquard
    %-'fmin': Using the fminsearch
    %algorithm
    switch cnf_p_ACDC.fitting_fun_type
        case 'lsq'
            %trust-region-reflective
            %levenberg-marquardt
            %'MaxIter',400
            cnf_p_ACDC.lsq_algorithm      = 'trust-region-reflective';
            cnf_p_ACDC.fitting_options    = optimset('Algorithm',cnf_p_ACDC.lsq_algorithm,'FinDiffType','central','TolFun',1e-12,'TolX',1e-12,'Display','off','UseParallel','always');
            %cnf_p_ACDC.fitting_options    = optimset('Algorithm',cnf_p_ACDC.lsq_algorithm,'Display','off','UseParallel','always');
            cnf_p_ACDC.fitting_options_lb = [0,0,0]; %lower boundaries
            cnf_p_ACDC.fitting_options_ACDC_lb = [-50.0,0,0];
            cnf_p_ACDC.fitting_options_ub = []; %upper boundaries
        case 'fmin'
            cnf_p_ACDC.fitting_options     =   optimset('Display','off','TolFun',1e-7,'TolX',1e-7);
            cnf_p_ACDC.fitting_options_lb = []; %lower boundaries
            cnf_p_ACDC.fitting_options_ACDC_lb = [-50.0,0,0];
            cnf_p_ACDC.fitting_options_ub = []; %upper boundaries
    end
    
    cnf_p_ACDC.geo_corr_application_flag=1;
    cnf_p_ACDC.force_geocorr_surf_type=1; % 0: use the geophysical corrections depending on the flag surface type included in the L1B product
    % 1: use the same group of geophysical corrections to all the records independently of the flag of the surface type: forcing the type surface
    cnf_p_ACDC.product_type_surface='open_ocean'; % open_ocean, sea_ice, land_ice ...
    
    %--------------- Verbose ploting --------------------------------------
    cnf_p_ACDC.plot_fits_flag    =   1;              % 0 - not plot; 1 - plot fit results
    cnf_p_ACDC.plot_fits_lat_range = [-91,91]; %latitude range of waveforms fit to be plotted (by default between -91 and 91 to force all of them to be plotted)
    cnf_p_ACDC.plot_fits_downsampling =100;
end
