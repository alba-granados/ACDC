% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the Instrument Gain Correction as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: The purpose of the instrument gain corrections is to compute 
% the gain/attenuation applied to the SAR science waveforms
% 
% INPUTs: ISP.power_val_cal1_sar_ku
%  -ISP.power_val_cal1_sar_ku                      dB	do  ISP decoding (§3.1)
%  -agc_telem_to_meas_ table_agc_cal    dB	do  AGC CAL selection (§5.1)
%  -power_cor_cal1                      dB	do  CAL1 selection (§5.1)
%  -ins_losses_ground_chd               dB	do  CHD file (§2.2.1)
%  -onboard_proc_sar_chd                dB	do  CHD file (§2.2.1)
%  -tx_rx_gain_ground_chd               dB	do  CHD file (§2.2.1)
%          
% OUTPUTs:  
%  -L1A.att_conv_sar_ku                     dB	do	
%  -L1A.gain_corr_instr_sar                 dB	do	
%  -power_var_cal1_sar_ku               dB	do	
%
% ----------------------------------------------------------
% Author:    Roger Escolà  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Mònica Roca   / isardSAT
% Last rev.: Mònica Roca   / isardSAT (11/09/2013)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L1A] = instrument_gain_corr (ISP)
global ins_losses_ground_chd power_var_cal1_sar_ku_rep onboard_proc_sar_chd rfu_rx_gain_ground_chd 

L1A.att_conv_sar_ku         = zeros(1,ISP.N_total_bursts_sar_ku);
L1A.gain_corr_instr_sar     = zeros(1,ISP.N_total_bursts_sar_ku);
L1A.power_val_cal1_sar_ku   = zeros(1,ISP.N_total_bursts_sar_ku);

for i_burst = 1:ISP.N_total_bursts_sar_ku

    progressbar([],[],[],[],[],(i_burst)/ISP.N_total_bursts_sar_ku,[],[],[],[],[],[],[],[],[]);

    % AGC tables adressing
    % L1A.att_conv_sar_ku = agc_telem_to_meas_table_cal1_chd (ISP.power_val_cal1_sar_ku);
    L1A.att_conv_sar_ku(i_burst) = ISP.att_sar_ku_isp(i_burst);

    % Instrument corrections
    L1A.gain_corr_instr_sar(i_burst) = rfu_rx_gain_ground_chd - L1A.att_conv_sar_ku(i_burst) - ins_losses_ground_chd + power_var_cal1_sar_ku_rep + onboard_proc_sar_chd;

    L1A.power_val_cal1_sar_ku(i_burst) = power_var_cal1_sar_ku_rep;
    
end
end
