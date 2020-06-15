%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the read RMC matrix
% algorithm as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: Read RMC data
% 
% INPUTs : Filename to read
% OUTPUTs: RMC matrix directly to be applied
%
% ----------------------------------------------------------
% Author:    Albert Garcia / isardSAT
%            Roger Escolà  / isardSAT
% Reviewer:  Mònica Roca   / isardSAT
% Last rev.: Mònica Roca   / isardSAT (11/09/2013)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [rmc_matrix] = read_RMC_matrix(filename_RMC)

global N_ku_pulses_burst_chd
global N_samples_sar_chd

fid_RMC  = fopen(filename_RMC,'r','b');
rmc_matrix_I = zeros(N_samples_sar_chd,N_ku_pulses_burst_chd);
rmc_matrix_Q = zeros(N_samples_sar_chd,N_ku_pulses_burst_chd);

for i_pulse = 1: N_ku_pulses_burst_chd 
    for i_sample = 1:N_samples_sar_chd
        rmc_matrix_I(i_sample,i_pulse) = 10^-4*fread(fid_RMC,1,'int16');% 256x64 [samples x pulses] elements with each element composed by I&Q samples stored as signed integers in 16+16 bits 
     end
end
for i_pulse = 1: N_ku_pulses_burst_chd 
    for i_sample = 1:N_samples_sar_chd
        rmc_matrix_Q(i_sample,i_pulse) = 10^-4*fread(fid_RMC,1,'int16');% 256x64 [samples x pulses] elements with each element composed by I&Q samples stored as signed integers in 16+16 bits 
    end
end
rmc_matrix= ((rmc_matrix_I + 1i*rmc_matrix_Q));
% figure;imagesc((abs((rmc_matrix))).');
% figure;imagesc((abs((fft(rmc_matrix)))/N_samples_sar_chd).');
fclose(fid_RMC);
end