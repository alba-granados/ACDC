function [burst] = coregistration_channels_SARin(burst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% DeDop Matlab Code
% 
% This code implements the algorithm to perform co-registration or
% alignment of the two burst of the two interferometric channels based on
% the window delay information (common window delay): coregistration done to reference channel-1
% used as reference during the whole processing: surface locations,...
%
% ---------------------------------------------------------
%
% Calling
%   burst = coregistration_channels_SARin( burst)
%
% Inputs
%   burst: input burst structure
%
% Output
%   modified burst          
% ----------------------------------------------------------
% 
% Author:   Eduard Makhoul / isardSAT
%
% Version  record
% 1.0 2016/09/22 Based on the window delay difference
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N_samples N_ku_pulses_burst_chd pi_cst

i_samples               = (0:(N_samples-1));

%% ------------------------------ Alignment --------------------
% Alignment using the phase ramp information similar to geometric
% corrections
%shift based on window delay difference
shift=(burst.win_delay_sar_ku-burst.win_delay_sar_ku_2)/burst.T0_sar;
burst.wfm_cal_gain_corrected_2  = burst.wfm_cal_gain_corrected_2.* ...
                    exp(-1i*2*pi_cst*...
                    (shift.*ones(N_ku_pulses_burst_chd,N_samples)).*...
                    (ones(N_ku_pulses_burst_chd,1)*i_samples/N_samples)); 

end