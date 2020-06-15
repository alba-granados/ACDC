%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% FAI APPLICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this function is to perform the alignment with the 
% FAI information assuming it is not performe on-board.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%____________________________For each stack________________________________

function [burst] = FAI_alignment (burst,FAI)

global N_samples N_ku_pulses_burst_chd c_cst pi_cst
global h0_cor2_unit_conv_chd T0_h0_unit_conv_chd

i_samples               = (0:(N_samples-1));

%% ------------------------------ FAI alignment ---------------------------
% In CryoSAt-2 no FAI alignment has been performed
% Compensate only for CAI
burst.wfm_cal_gain_corrected  = burst.wfm_cal_gain_corrected.* ...
                    exp(1i*2*pi_cst*...
                    (((FAI-FAI(1))./ h0_cor2_unit_conv_chd / T0_h0_unit_conv_chd).'*ones(1,N_samples)).*...
                    (ones(N_ku_pulses_burst_chd,1)*i_samples/N_samples));

% Need to include the second channel for SARIn
       

end