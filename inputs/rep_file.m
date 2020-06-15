% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% CR2 SARIn 
%
% ---------------------------------------------------------
% Objective: Define Repositori Parameters
% 
% INPUTs : - 
% OUTPUTs: -
%
% ----------------------------------------------------------
% Author:    Roger Escolà  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Mònica Roca   / isardSAT
% Last rev.: Mònica Roca   / isardSAT (11/09/2013)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global int_delay_cal1_sar_iq_ku_rep
global power_var_cal1_sar_ku_rep burst_phase_array_cor_cal1_sar_rep
global burst_power_array_cor_cal1_sar_rep wfm_cal2_science_sar_rep 
global N_ku_pulses_sar_chd N_samples_sar_chd


int_delay_cal1_sar_iq_ku_rep = 10 * 1e-9;
power_var_cal1_sar_ku_rep = 0;
burst_phase_array_cor_cal1_sar_rep = zeros(1,N_ku_pulses_sar_chd);
burst_power_array_cor_cal1_sar_rep = zeros(1,N_ku_pulses_sar_chd);
wfm_cal2_science_sar_rep = ones (1,N_samples_sar_chd);
