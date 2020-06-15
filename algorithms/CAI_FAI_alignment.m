function [burst] = CAI_FAI_alignment (burst,CAI,FAI)

global N_samples N_ku_pulses_burst_chd c_cst pi_cst
global h0_cor2_unit_conv_chd T0_h0_unit_conv_chd cai_cor2_unit_conv_chd
global mode

i_samples               = (0:(N_samples-1));

%% ------------------------------ FAI alignment ---------------------------
% In CryoSAt-2 no FAI alignment has been performed
% Compensate only for CAI

%shift related to the difference in window delay due to internal path
%delay corrections
burst.wfm_cal_gain_corrected  = burst.wfm_cal_gain_corrected.* ...
                    exp(-1i*2*pi_cst*...
                    (((FAI)./ h0_cor2_unit_conv_chd / T0_h0_unit_conv_chd).'*ones(1,N_samples)).*...
                    (ones(N_ku_pulses_burst_chd,1)*i_samples/N_samples));
% burst.wfm_cal_gain_corrected_FAI_positive_sign  = burst.wfm_cal_gain_corrected_before_align.* ...
%                     exp(1i*2*pi_cst*...
%                     ((((CAI-CAI(1)).*cai_cor2_unit_conv_chd+FAI)./ h0_cor2_unit_conv_chd / T0_h0_unit_conv_chd).'*ones(1,N_samples)).*...
%                     (ones(N_ku_pulses_burst_chd,1)*i_samples/N_samples));                

if strcmp(mode,'SIN')
   burst.wfm_cal_gain_corrected_2  = burst.wfm_cal_gain_corrected_2.* ...
                    exp(-1i*2*pi_cst*...
                    (((FAI)./ h0_cor2_unit_conv_chd / T0_h0_unit_conv_chd).'*ones(1,N_samples)).*...
                    (ones(N_ku_pulses_burst_chd,1)*i_samples/N_samples)); 
end
%include the impact of FAI in the window delay
burst.win_delay_sar_ku=burst.win_delay_sar_ku+mean(((CAI-CAI(1)).*cai_cor2_unit_conv_chd+(FAI))./ h0_cor2_unit_conv_chd / T0_h0_unit_conv_chd.*burst.T0_sar);       


% %% Validation purposes
% % Detection of leading edge
% percent_leading_edge=75.0;
% wfm_before_az_noFAI_align=abs(fftshift(fft(burst.wfm_cal_gain_corrected_before_align,[],2))).^2;
% wfm_before_az_FAI_align=abs(fftshift(fft(burst.wfm_cal_gain_corrected,[],2))).^2;
% wfm_before_az_FAI_align_plusign=abs(fftshift(fft(burst.wfm_cal_gain_corrected_FAI_positive_sign,[],2))).^2;
% [peak_pow_noFAI,idx_max_peak_noFAI]=max(wfm_before_az_noFAI_align,[],2);
% [peak_pow_FAI,idx_max_peak_FAI]=max(wfm_before_az_FAI_align,[],2);
% [peak_pow_FAI_plussign,idx_max_peak_FAI_plussign]=max(wfm_before_az_FAI_align_plusign,[],2);
% for i_wfm=1:N_ku_pulses_burst_chd
%     %before FAI
%     dumm=find(find(wfm_before_az_noFAI_align(i_wfm,:)<=percent_leading_edge/100.0*peak_pow_noFAI(i_wfm))<idx_max_peak_noFAI(i_wfm), 1, 'last' );
%     if ~isempty(dumm)
%         burst.idx_leading_noFAI(i_wfm)=dumm;
%     else
%         %if there is no leading edge or the waveform has displaced that
%         %much from the window to the left select the peak as leading
%         %edge
%         burst.idx_leading_noFAI(i_wfm)=idx_max_peak_noFAI(i_wfm);
%     end
%     
%     %after FAI
%     dumm=find(find(wfm_before_az_FAI_align(i_wfm,:)<=percent_leading_edge/100.0*peak_pow_FAI(i_wfm))<idx_max_peak_FAI(i_wfm), 1, 'last' );
%     if ~isempty(dumm)
%         burst.idx_leading_FAI(i_wfm)=dumm;
%     else
%         %if there is no leading edge or the waveform has displaced that
%         %much from the window to the left select the peak as leading
%         %edge
%         burst.idx_leading_FAI(i_wfm)=idx_max_peak_FAI(i_wfm);
%     end
%     
%     %after FAI PLUS SIGN
%     dumm=find(find(wfm_before_az_FAI_align_plusign(i_wfm,:)<=percent_leading_edge/100.0*peak_pow_FAI_plussign(i_wfm))<idx_max_peak_FAI_plussign(i_wfm), 1, 'last' );
%     if ~isempty(dumm)
%         burst.idx_leading_FAI_plussign(i_wfm)=dumm;
%     else
%         %if there is no leading edge or the waveform has displaced that
%         %much from the window to the left select the peak as leading
%         %edge
%         burst.idx_leading_FAI_plussign(i_wfm)=idx_max_peak_FAI_plussign(i_wfm);
%     end
% 
% end
end