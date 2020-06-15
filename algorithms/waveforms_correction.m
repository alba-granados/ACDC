% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the WAVEFORMS CORRECTION as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: - Correct every echo of a burst from power and phase variations 
% within a burst using instrumental calibration corrections (CAL1).
%  -Correct waveforms for the CAL2 provided by instrumental calibration
%   corrections (CAL2).
% 
% INPUTs: 
%   
%  -N_ku_pulses_burst_chd                         uc	CHD file (§2.2.1)
%  -N_samples_sar_chd                           us	CHD file (§2.2.1)
%  -fai_fine_shift_number_chd                   us	CHD file (§2.2.1)
%  -L1A.gain_corr_instr_sar                     dB	do	Instrument Gain (§7.4)
%  -burst_phase_array_cor_cal1_sar          rad	do	CAL1 selection (§5.1)
%  -burst_power_array_cor_cal1_sar_rep      dB	do	CAL1 selection (§5.1)
%  -wfm_cal2_science_sar                    FFT ss	CAL2 selection (§5.2)
%  -fai_sar                                 s	do	Preliminary Window Delay (§7.3)
%  -ISP.wfm_sar_reversed                        FFT do  OnBoard reversion (§7.43.1)
%          
% OUTPUTs:  
%  -burst_phase_array_cor_cal1_sar          rad do  
%  -burst_power_array_cor_cal1_sar          dB  do
%  -wfm_cal2_science_sar                    FFT do 
%  -wfm_cal_corrected                       FFT do 
%
% ----------------------------------------------------------
% Author:    Roger Escolà  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Mònica Roca   / isardSAT
% Last rev.: Mònica Roca   / isardSAT (11/09/2013)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L1A] = waveforms_correction (ISP, L1A)

global burst_power_array_cor_cal1_sar_rep burst_phase_array_cor_cal1_sar_rep wfm_cal2_science_sar_rep 
global N_ku_pulses_burst_chd N_samples_sar_chd                                                


wfm_sar_ku_cal1         = zeros(ISP.N_total_bursts_sar_ku,N_ku_pulses_burst_chd,N_samples_sar_chd);
wfm_gain_corrected      = zeros(ISP.N_total_bursts_sar_ku,N_ku_pulses_burst_chd,N_samples_sar_chd);
wfm_cal_corrected_aux   = zeros(N_samples_sar_chd,N_ku_pulses_burst_chd);
L1A.wfm_cal_gain_corrected  = zeros(ISP.N_total_bursts_sar_ku,N_ku_pulses_burst_chd,N_samples_sar_chd);

for i_burst = 1:ISP.N_total_bursts_sar_ku
    
    progressbar([],[],[],[],[],[],(i_burst-1)/ISP.N_total_bursts_sar_ku,[],[],[],[],[],[],[],[]);
   
    % Putting together the real and the imaginary parts
    % aux (1:N_ku_pulses_burst_chd, 1:N_samples_sar_chd) = wfm_iq_sar_ku_isp (1,:,:) + 1i * wfm_iq_sar_ku_isp (2,:,:);
    % ALREADY DONE in the reading of the FBR!!

    
    wfm_gain_corrected(i_burst,:,:) = ISP.wfm_sar_reversed(i_burst,:,:) ./ 10^(L1A.gain_corr_instr_sar(i_burst)/20);
    
    % CAL1 Application
    for i_pulse = 1:N_ku_pulses_burst_chd
        wfm_sar_ku_cal1 (i_burst,i_pulse,:) = wfm_gain_corrected (i_burst,i_pulse,:) * (10 ^ (burst_power_array_cor_cal1_sar_rep(i_pulse)/20)) * exp (1i * burst_phase_array_cor_cal1_sar_rep(i_pulse));
        % the way burst_phase_array_cor_cal1_sar_rep has been computed has to
        % be changed to " - [phase(i) - phase(1)] "
    end

    %***************************************************
    % CAL2 Application
    % wfm_cal2_science_sar_rep MUST be a column vector

    % Swap the CAL2 waveform if necessary %
    % SWAP:
    % cal2_wfm = fftshift (wfm_cal2_science_sar_rep);

    % NO SWAP:
    cal2_wfm = wfm_cal2_science_sar_rep;
    %***************************************************


    % FFT
    aux_burst(1:N_ku_pulses_burst_chd,1:N_samples_sar_chd) = wfm_sar_ku_cal1(i_burst,:,:);
    wfm_fft_sar_ku_cal1 = fft(aux_burst.'); % The result is in columns

    for i_pulse = 1:N_ku_pulses_burst_chd
        % The multiplication is performed in columns (because of the FFT's result)
        wfm_cal_corrected_aux (:, i_pulse) = wfm_fft_sar_ku_cal1 (:, i_pulse) .* cal2_wfm.';
    end

    % IFFT
    L1A.wfm_cal_gain_corrected(i_burst,:,:) = ifft(wfm_cal_corrected_aux).';

end
end
