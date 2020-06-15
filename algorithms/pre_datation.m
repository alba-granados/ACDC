% PRELIMINARY DATATION ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% DeDop 
% This code implements the PRELIMINARY DATATION as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5b_20131125
%
% ---------------------------------------------------------
% Objective: Compute the preliminary waveform % datation
% 
% ----------------------------------------------------------
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/09/2013)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L1A]   = pre_datation (ISP)

%% 0.Initialisation and Memory allocation

global pulse_length_chd tx1_sar_chd pri_T0_unit_conv_chd T0_nom
global mean_sat_alt_cst c_cst
global N_total_bursts_sar_ku_isp N_pulses_burst
global bri_nom i_ISP



tp                          = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
t_tx                        = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));
pri_sar_nom                 = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));
L1A.time_sar_ku_pre_dat     = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));
L1A.time_sar_ku_pre_dat_pulse   = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
burst_prop_delay            = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));

for i_burst = 1:N_total_bursts_sar_ku_isp(i_ISP)
    progressbar([],(i_burst)/N_total_bursts_sar_ku_isp(i_ISP),[],[],[],[],[],[],[],[],[],[],[],[]);


    %% 1.Time stamp delay
    burst_prop_delay(i_burst) = bri_nom(i_burst) * (ISP.burst_sar_isp(i_burst)-1); %this is the propagation of the time stamp along the radar cycle
    t_tx(i_burst) = ISP.time_sar_isp(i_burst) + tx1_sar_chd + pulse_length_chd/2 + burst_prop_delay(i_burst);

    %% 2.Shift to each transmitted pulse
    pri_sar_nom(i_burst) = ISP.pri_sar_isp(i_burst) * pri_T0_unit_conv_chd * T0_nom;
    
    for i_pulse = 1:N_pulses_burst
        tp(i_burst,i_pulse) = t_tx(i_burst) - (ISP.ambiguity_order_sar_isp(i_burst)-i_pulse+1) * pri_sar_nom(i_burst);
    end

    %% 3.Locate at the surface
    delta_prop_delay                        = mean_sat_alt_cst/c_cst;
    L1A.time_sar_ku_pre_dat_pulse(i_burst,:)    = tp(i_burst,:) + delta_prop_delay;

    %% 4.Datation averaging
    L1A.time_sar_ku_pre_dat(i_burst) = mean(L1A.time_sar_ku_pre_dat_pulse(i_burst,2:N_pulses_burst));
    
end


end