%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% HEIGHT RATE APPLICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this function is to perform a better alignment with the 
% height rate information once the applied CAI/FAI have been compensated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v1.1 Changed N_samples > N_samples

%____________________________For each stack________________________________

function [burst] = height_rate_alignment (burst,CAI)

global N_samples N_ku_pulses_burst_chd c_cst pi_cst
global cai_cor2_unit_conv_chd h0_cor2_unit_conv_chd T0_h0_unit_conv_chd
global mode

i_samples               = (0:(N_samples-1));
i_pulses                = (0:(N_ku_pulses_burst_chd-1)).'; 

% %% ------------------------------ CAI/FAI compensation --------------------
% % In CryoSAt-2 no FAI alignment has been performed
% % Compensate only for CAI
% burst.wfm_cal_gain_corrected  = burst.wfm_cal_gain_corrected.* ...
%                     exp(-1i*2*pi_cst*...
%                     (((CAI-CAI(1)).*cai_cor2_unit_conv_chd / h0_cor2_unit_conv_chd / T0_h0_unit_conv_chd)*ones(1,N_samples)).*...
%                     (ones(N_ku_pulses_burst_chd,1)*i_samples/N_samples));

%% ----------------------------- Height rate application ------------------
burst.wfm_cal_gain_corrected  = burst.wfm_cal_gain_corrected.* ...
                    exp(-1i*2*pi_cst*burst.alt_rate_sar_sat.*...
                    ((i_pulses)*ones(1,N_samples)).*burst.pri_sar * 2/c_cst/burst.T0_sar.*...
                    (ones(N_ku_pulses_burst_chd,1)*i_samples/N_samples));
                
if strcmp(mode,'SIN')
    burst.wfm_cal_gain_corrected_2  = burst.wfm_cal_gain_corrected.* ...
                    exp(-1i*2*pi_cst*burst.alt_rate_sar_sat.*...
                    ((i_pulses)*ones(1,N_samples)).*burst.pri_sar * 2/c_cst/burst.T0_sar.*...
                    (ones(N_ku_pulses_burst_chd,1)*i_samples/N_samples));
    
end
%% ----------------------------- Window delay update ----------------------
alt_rate_wd_corr = burst.alt_rate_sar_sat.*...
            (i_pulses) * burst.pri_sar * 2/c_cst; % in seconds
burst.win_delay_sar_ku=burst.win_delay_sar_ku+mean(alt_rate_wd_corr);        

end