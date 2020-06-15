function wfm_gain_CAL1_CAL2_corrected = waveforms_correction_S3(wfm_cal_gain_UNCORRECTED,netCDF_L1A)






% This algorithm follow DPM S3IPF.DPM.005 
% ALT_COR_WAV_05 AGC
% ALT_COR_WAV_02 Phase intra Burst
% ALT_COR_WAV_03 Power Intra Burst
% ALT_COR_WAV_04 CAL2
% The sig0_cal_ku_l1a_echo_sar_ku has been included although it is not applied anywhere in the S3DPM.

% It is called after the L1A record is read, during the adaptation fuction.

% Author:   Albert Garcia-Mondejar / isardSAT
% v1.0 first version of the algorithm. 



global N_samples_sar_chd N_ku_pulses_burst_chd
fixed_gain =    - double(netCDF_L1A.data.sig0_cal_ku_l1a_echo_sar_ku) .*double(netCDF_L1A.attributes.sig0_cal_ku_l1a_echo_sar_ku.scale_factor)+ ... % dB 'internal calibration correction on Sigma0 for ku band: l1a_echo_sar_ku mode'
                double(netCDF_L1A.data.agc_ku_l1a_echo_sar_ku)      .*double(netCDF_L1A.attributes.agc_ku_l1a_echo_sar_ku.scale_factor); % dB corrected AGC for ku band: l1a_echo_sar_ku mode


%% CAL1 Application
CAL1_power                  = double(netCDF_L1A.data.burst_power_cor_ku_l1a_echo_sar_ku)    .*double(netCDF_L1A.attributes.burst_power_cor_ku_l1a_echo_sar_ku.scale_factor) *ones(1,N_samples_sar_chd); %'ku band burst power corrections (cal1) : l1a_echo_sar_ku mode' %'FFT power unit'
CAL1_phase                  = double(netCDF_L1A.data.burst_phase_cor_ku_l1a_echo_sar_ku)    .*double(netCDF_L1A.attributes.burst_phase_cor_ku_l1a_echo_sar_ku.scale_factor) *ones(1,N_samples_sar_chd);%'ku band burst phase corrections (cal1) : l1a_echo_sar_ku mode' %radians
CAL2                        = double(netCDF_L1A.data.gprw_meas_ku_l1a_echo_sar_ku(:,3))     .*double(netCDF_L1A.attributes.gprw_meas_ku_l1a_echo_sar_ku.scale_factor)       *ones(1,N_ku_pulses_burst_chd); %128x3 OPTION 3 chosen as it seems less noisier.  Normalized GPRW (cal2) %'FFT power unit'

wfm_gain_CAL1_corrected     = wfm_cal_gain_UNCORRECTED .* (sqrt(CAL1_power).*exp (1i .* CAL1_phase)) ./ sqrt(10.^(fixed_gain/10));



%% CAL2 Application


wfm_gain_CAL1_corrected_fft = fftshift(fft(wfm_gain_CAL1_corrected.'));
wfm_gain_CAL1_CAL2_corrected_fft = wfm_gain_CAL1_corrected_fft ./ sqrt(CAL2);
wfm_gain_CAL1_CAL2_corrected = ifft(fftshift(wfm_gain_CAL1_CAL2_corrected_fft)).';

% figure; subplot(2,2,1); plot(sum(abs(fftshift(fft(wfm_gain_CAL1_CAL2_corrected.'))).'));
% hold all; plot(sum(abs(fftshift(fft(wfm_gain_CAL1_corrected.'))).'));
% hold all; plot(sum(abs(fftshift(fft(wfm_cal_gain_UNCORRECTED.'))).'));
% legend('CAL1 and CAL2  corrected','CAL1   corrected', 'L1A waveforms uncorrected');
% 
% subplot(2,2,2); plot(double(netCDF_L1A.data.burst_power_cor_ku_l1a_echo_sar_ku)    .*double(netCDF_L1A.attributes.burst_power_cor_ku_l1a_echo_sar_ku.scale_factor));
% figlabels('Pulse Index','FFT units','','CAL1 power',12);
% subplot(2,2,4); plot(double(netCDF_L1A.data.burst_phase_cor_ku_l1a_echo_sar_ku)    .*double(netCDF_L1A.attributes.burst_phase_cor_ku_l1a_echo_sar_ku.scale_factor));
% figlabels('Pulse Index','radians','','CAL1 phase',12);
% subplot(2,2,3); plot(double(netCDF_L1A.data.gprw_meas_ku_l1a_echo_sar_ku(:,3))     .*double(netCDF_L1A.attributes.gprw_meas_ku_l1a_echo_sar_ku.scale_factor));
% figlabels('Range bin','FFT units','','CAL2',12);



end