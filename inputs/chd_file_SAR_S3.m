% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% CR2 SARIn %
% ---------------------------------------------------------
% Objective: Define characterization parameters
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Other (not in the DPM)
% global alt_clock_period_sar_ku
global c_cst 
global date_ref_switch_SIRAL date_ref_time_SIRAL


% alt_clock_period_sar_ku = 1.25 * 1e-8; % 1/f0


date_ref_switch_SIRAL=datevec('October 21, 2010 00:00:00','mmmm dd, yyyy HH:MM:SS'); %assume is 00:00:00 hours
date_ref_time_SIRAL=datevec('January 1, 2000 00:00:00','mmmm dd, yyyy HH:MM:SS');
date_ref_switch_SIRAL=etime(date_ref_switch_SIRAL,date_ref_time_SIRAL); %w.r.t January January 1, 2000 00:00:00 in seconds


%% Main
global freq_ku_chd bw_ku_chd wv_length_ku
global chirp_slope_ku_chd

freq_ku_chd             = 13.575 * 1e9; % c_cst / wv_length_ku;
bw_ku_chd               = 320 * 1e6;
bw_c_chd                = 290 * 1e6;
wv_length_ku 			= c_cst / freq_ku_chd;%
% freq_c_chd              = c_cst / wv_length_c;
% bw_c_chd                = 'TBC' * 1e6;

%% Time patern
global N_samples_sar_chd
global N_ku_pulses_sar_chd  N_c_pulses_burst_chd 
global N_bursts_cycle_chd N_pri_sar_c_ku_chd
global tx1_sar_chd pulse_length_chd burst_duration_sar_chd 
global prf_sar_chd brf_chd N_ku_pulses_burst_chd
global N_max_beams_stack_chd pri_sar_chd bri_chd
global PTR_width_chd

N_samples_sar_chd       = 128;    % SAR Samples

N_ku_pulses_burst_chd   = 64;   % SAR Ku pulses in burst
N_c_pulses_burst_chd    = 1;     % SAR C pulses in burst
N_ku_pulses_sar_chd     = N_ku_pulses_burst_chd + N_c_pulses_burst_chd;
N_bursts_cycle_chd  = 4; % Bursts in a cycle
N_pri_sar_c_ku_chd      = 0; %NOTFOUND
N_max_beams_stack_chd   = N_bursts_cycle_chd*N_ku_pulses_burst_chd;

pulse_length_chd        = 44.799999 * 1e-6;
tx1_sar_chd             = 0; %NOTFOUND -3.7116047 * 1e-3; %cog_antenna_dat_error

prf_sar_chd = 17.825311 * 1e+03; %[Hz]
brf_chd = 78.5306931;

pri_sar_chd = 1 / prf_sar_chd 
bri_chd = 1 / brf_chd; %11.7929625 * 1e-3*4;

chirp_slope_ku_chd      = bw_ku_chd/pulse_length_chd;
burst_duration_sar_chd  = N_ku_pulses_sar_chd * pri_sar_chd;
% burst_duration_sar_chd = 0.0035904;


%to be used in the sigma0 to be in accordance with L2ESA/technical note
%3dBbeamwidth
PTR_width_chd = 19821; %NOTSURE %2.819e-9; %according to Dinardo's note 2.819e-9 seconds instead of the 1/(chirp_slope_ku_chd*BW) 
%% Platform
global x_ant_chd y_ant_chd z_ant_chd x_cog_chd y_cog_chd z_cog_chd
global x_cog_ant y_cog_ant z_cog_ant

x_ant_chd               = 0;
y_ant_chd               = 0;
z_ant_chd               = 0;
x_cog_chd               = 0;
y_cog_chd               = 0;
z_cog_chd               = 0;

x_cog_ant = x_ant_chd - x_cog_chd;
y_cog_ant = y_ant_chd - y_cog_chd;
z_cog_ant = z_ant_chd - z_cog_chd;

%% Antenna
global antenna_gain_ku_chd antenna_beamwidth_ku_chd
%global roll_random_error_chd pitch_random_error_chd yaw_random_error_chd
%global roll_bias_error_chd pitch_bias_error_chd yaw_bias_error_chd
%global roll_harmonic_error_chd pitch_harmonic_error_chd yaw_harmonic_error_chd
%Errors of pointing from Rob presentation on requirements meeting 28/01/2015

%roll_random_error_chd       = 0.0496/180*pi; %radians.
%pitch_random_error_chd      = 0.0499/180*pi; %radians.
%yaw_random_error_chd        = 0.0494/180*pi; %radians.
%roll_bias_error_chd         = 0.0828/180*pi; %radians.
%pitch_bias_error_chd        = 0.0722/180*pi; %radians.
%yaw_bias_error_chd          = 0.0068/180*pi; %radians.
%roll_harmonic_error_chd     = 0.0441/180*pi; %radians.
%pitch_harmonic_error_chd    = 0.0534/180*pi; %radians.
%yaw_harmonic_error_chd      = 0.0250/180*pi; %radians.

antenna_gain_ku_chd             = 41.9000015;
antenna_beamwidth_ku_chd        = 1.35; %Degrees

%% Window Delay
global ext_delay_ground_chd 
%global int_delay_ground_chd


ext_delay_ground_chd = 1.1 * 1e-05; %15386 * 1e-12;
%int_delay_ground_chd = - 7000e-12;


%% AGC and waveforms scaling factor computation
global agc_telem_to_meas_table_cal1_chd 
global onboard_proc_sar_chd 
global ins_losses_ground_chd power_tx_ant_ku_chd
global rfu_rx_gain_ground_chd
global ADC_mult_factor PTR_power_drift_slope_sec

agc_telem_to_meas_table_cal1_chd = 0; %NOTFOUND
rfu_rx_gain_ground_chd = 99; %NOTFOUND
ADC_mult_factor=1000; %NOTFOUND
onboard_proc_sar_chd = 10*log10(ADC_mult_factor.^2); % Assuming FFT are unitary transformations                                                                                                              % Compensate for pulse compression in range and DOppler (I don't agree) as Dinardo says %NOTFOUND
%onboard_proc_sar_chd =10*log10(ADC_mult_factor.^2)+10*log10(bw_ku_chd*pulse_length_chd); % Processing gain only on pulse compression: TBP
PTR_power_drift_slope_sec=-0.022/(30*24*60*60); %from Dinardo note


ins_losses_ground_chd = 0; %NOTFOUND
power_tx_ant_ku_chd = 20; %NOTFOUND

%% USO CLock
global uso_freq_nom_chd alt_freq_multiplier_chd T0_chd
%global T0_nom i_sample_start_chd
%global pri_T0_unit_conv_chd h0_cor2_unit_conv_chd cai_cor2_unit_conv_chd T0_h0_unit_conv_chd


uso_freq_nom_chd = 10e6;
alt_freq_multiplier_chd = 32; %NOTFOUND
T0_chd = 1/(uso_freq_nom_chd * alt_freq_multiplier_chd);

%pri_T0_unit_conv_chd = 8;
%h0_cor2_unit_conv_chd = 16;
%T0_h0_unit_conv_chd = 64;
%cai_cor2_unit_conv_chd = 4096;

%i_sample_start_chd = 1;




%% Weightings
%azimuth_weighting_filename_chd = [];











