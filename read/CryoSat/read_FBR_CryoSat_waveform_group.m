function [burst] = read_FBR_CryoSat_waveform_group(fid, burst)

global mode N_samples
global N_ku_pulses_burst_chd N_samples_sar_chd N_samples_sin_chd    
global burst_phase_array_cor_cal1_sar_rep_1 burst_phase_array_cor_cal1_sar_rep_2 burst_power_array_cor_cal1_sar_rep_1 burst_power_array_cor_cal1_sar_rep_2 
global wfm_cal2_science_sar_rep_1 wfm_cal2_science_sar_rep_2
global CAL1p2p_flag_cnf CAL2_flag_cnf
global height_rate_application_cnf FAI_application_cnf

%% INITIALISATION
switch mode
    case 'SAR'
        N_samples = N_samples_sar_chd;
        
    case 'SIN'
        N_samples = N_samples_sin_chd;
        
        wfm_iq_sar_ku_fbr_aux_2 = zeros(1,N_ku_pulses_burst_chd*N_samples*2);
        wfm_cal_gain_corrected_2_i  = zeros(N_ku_pulses_burst_chd,N_samples);
        wfm_cal_gain_corrected_2_q  = zeros(N_ku_pulses_burst_chd,N_samples);
end

wfm_iq_sar_ku_fbr_aux_1 = zeros(1,N_ku_pulses_burst_chd*N_samples*2);
wfm_cal_gain_corrected_1_i  = zeros(N_ku_pulses_burst_chd,N_samples);
wfm_cal_gain_corrected_1_q  = zeros(N_ku_pulses_burst_chd,N_samples);



%   fread(fid,2621520,'int8'); 
%   SAR FBR waveform group 327760  
%   SARIn FBR waveforms group 2621520 
    
        
%% READING
wfm_iq_sar_ku_fbr_aux_1(:) = fread (fid,N_ku_pulses_burst_chd*N_samples*2,'int8');
if strcmp(mode,'SIN')
    wfm_iq_sar_ku_fbr_aux_2(:) = fread (fid,N_ku_pulses_burst_chd*N_samples*2,'int8');
end
wfm_iq_sar_ku_fbr_1_i = wfm_iq_sar_ku_fbr_aux_1(1:2:N_ku_pulses_burst_chd * N_samples*2);
wfm_iq_sar_ku_fbr_1_q = wfm_iq_sar_ku_fbr_aux_1(2:2:N_ku_pulses_burst_chd * N_samples*2);

if strcmp(mode,'SIN')
    wfm_iq_sar_ku_fbr_2_i = wfm_iq_sar_ku_fbr_aux_2(1:2:N_ku_pulses_burst_chd * N_samples*2);
    wfm_iq_sar_ku_fbr_2_q = wfm_iq_sar_ku_fbr_aux_2(2:2:N_ku_pulses_burst_chd * N_samples*2);
end

        
%% APPLICATION OF INST GAIN CORR
for i_pulse = 1:N_ku_pulses_burst_chd
    wfm_cal_gain_corrected_1_i(i_pulse,:) = wfm_iq_sar_ku_fbr_1_i((i_pulse-1)*N_samples+1:i_pulse*N_samples)./(10.^(burst.gain_corr_instr_sar_1./20));
    wfm_cal_gain_corrected_1_q(i_pulse,:) = wfm_iq_sar_ku_fbr_1_q((i_pulse-1)*N_samples+1:i_pulse*N_samples)./(10.^(burst.gain_corr_instr_sar_1./20));
end
if strcmp(mode,'SIN')
    for i_pulse = 1:N_ku_pulses_burst_chd
        wfm_cal_gain_corrected_2_i(i_pulse,:) = wfm_iq_sar_ku_fbr_2_i((i_pulse-1)*N_samples+1:i_pulse*N_samples)./(10.^(burst.gain_corr_instr_sar_2./20));
        wfm_cal_gain_corrected_2_q(i_pulse,:) = wfm_iq_sar_ku_fbr_2_q((i_pulse-1)*N_samples+1:i_pulse*N_samples)./(10.^(burst.gain_corr_instr_sar_2./20));
    end
end
clear wfm_iq_sar_ku_fbr_i_1 wfm_iq_sar_ku_fbr_q_1 wfm_iq_sar_ku_fbr_i_2 wfm_iq_sar_ku_fbr_q_2

burst.wfm_cal_gain_corrected        = wfm_cal_gain_corrected_1_i + 1i*wfm_cal_gain_corrected_1_q;
if strcmp(mode,'SIN')
    burst.wfm_cal_gain_corrected_2  = wfm_cal_gain_corrected_2_i + 1i*wfm_cal_gain_corrected_2_q;
end


%% APPLICATION OF CAL1 intra burst
if CAL1p2p_flag_cnf
    burst.wfm_cal_gain_corrected = burst.wfm_cal_gain_corrected.*((10.^ (burst_power_array_cor_cal1_sar_rep_1.'/20))*ones(1,N_samples)).*...
        exp (1i * (burst_phase_array_cor_cal1_sar_rep_1.'*ones(1,N_samples)));
    
    if strcmp(mode,'SIN')
        burst.wfm_cal_gain_corrected_2 = burst.wfm_cal_gain_corrected_2.*((10.^ (burst_power_array_cor_cal1_sar_rep_2.'/20))*ones(1,N_samples)).*...
        exp (1i * (burst_phase_array_cor_cal1_sar_rep_2.'*ones(1,N_samples)));
    end
end


%% APPLICATION OF CAL2
if CAL2_flag_cnf
    aux = fftshift(fft(burst.wfm_cal_gain_corrected,N_samples,2),2); % fft along range samples
    burst.wfm_cal_gain_corrected = ifft(fftshift(aux.*(ones(N_ku_pulses_burst_chd,1)*wfm_cal2_science_sar_rep_1),2),N_samples,2);
    clear aux;
    
    if strcmp(mode,'SIN')
        aux = fftshift(fft(burst.wfm_cal_gain_corrected_2,N_samples,2),2); % fft along range samples
        burst.wfm_cal_gain_corrected_2 = ifft(fftshift(aux.*(ones(N_ku_pulses_burst_chd,1)*wfm_cal2_science_sar_rep_2),2),N_samples,2);
        clear aux;
    end
end

%% --------------- Alignment pulses by FAI or height rate -----------------
if height_rate_application_cnf
    %Computation of the CAI
    [CAI,~]=CAI_FAI_retrieval(burst);
    %FAI is not applied on-board
    [burst] = height_rate_alignment(burst,CAI);
else
    if FAI_application_cnf
        %Computation of the FAI
        [CAI,FAI]=CAI_FAI_retrieval(burst);
        [burst]=CAI_FAI_alignment(burst,CAI,FAI);
    end
end


%% FURTHER READING
% 52 Echo Scale Factor (to scale echo to watts) - 4 sl
burst.nimp_sar_isp = fread(fid,1,'uint16');

fread(fid,1,'uint16');
    

    
% end
clear wfm_cal_gain_corrected_1_i wfm_cal_gain_corrected_1_q
clear wfm_cal_gain_corrected_2_i wfm_cal_gain_corrected_2_q


% figure; imagesc(abs(fftshift(fft(squeeze(wfm_iq_sar_ku_fbr_i_2(1,:,:)-1i*wfm_iq_sar_ku_fbr_q_2(1,:,:))'),1)))
% figure; imagesc(abs(fftshift(fft(fftshift(fft(squeeze(burst.wfm_cal_gain_corrected(:,:))))'),1))')


end