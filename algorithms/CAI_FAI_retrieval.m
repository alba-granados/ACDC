%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Computation of the CAI/FAI p2p %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this function is to extract the CAI/FAI components pulse 
% to pulse using the information of H0 and COR2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%____________________________For each stack________________________________

function [CAI,FAI] = CAI_FAI_retrieval (burst)

global N_ku_pulses_burst_chd h0_cor2_unit_conv_chd cai_cor2_unit_conv_chd T0_h0_unit_conv_chd
global bri_chd N_bursts_cycle_chd
%% -------------------------- Height in seconds ---------------------------
% assuming the PRI is the same for all the bursts

total_num_pulses_RC=round(N_bursts_cycle_chd*bri_chd/burst.pri_sar);
h=floor(burst.h0_comp_sar_isp*h0_cor2_unit_conv_chd+burst.cor2_comp_sar_isp.*(0:1:(N_ku_pulses_burst_chd-1))/(total_num_pulses_RC)); %units of COR2

%% ----------------------------- CAI extraction ---------------------------
CAI = floor(h./cai_cor2_unit_conv_chd); % units: T0*4

%% ----------------------------- FAI extraction ---------------------------
FAI = h - CAI.* cai_cor2_unit_conv_chd; % units: T0/64/16

end