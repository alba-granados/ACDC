function [burst] = read_FBR_CryoSat_measurements_group(fid, burst)

global power_tx_ant_ku_chd onboard_proc_sar_chd
%EM 13.04.2016
global zp_fact_range_cnf date_ref_switch_SIRAL PTR_power_drift_slope_sec
global gain_scale_hamm gain_scale_win_rg
global N_samples N_ku_pulses_burst_chd;
global bri_chd N_bursts_cycle_chd;
global mode
% global pulse_length_chd chirp_slope_ku_chd

    %-----------------------------%
    %-- read Measurements Group --%
    %-----------------------------%
                    
		%15 Window Delay (2way) uncorrected for instrument delays 10-12 s 8 sll
        burst.win_delay_sar_ku = fread(fid,1,'int64') * 1e-12;
		%16 Ho Initial Height Word 48.8 ps 4 sl (see note 2)
        burst.h0_comp_sar_isp  = fread(fid,1,'int32');% * 48.8 *1e-12;
        %17 HPR Height Rate 3.05 ps/rc 4 sl (see note 2)
        burst.cor2_comp_sar_isp = fread(fid,1,'int32');% * 3.05 * 1e-12/(bri_chd*N_bursts_cycle_chd);
		%18 LAI 12.5 ns 4 sl (see note 2)
        burst.cai_sar_isp = fread(fid,1,'int32');% * 12.5 * 1e-9;
		%19 FAI 12.5/256 ns 4 sl (see note 2)
        burst.fai_sar_isp = fread(fid,1,'int32');% * 12.5/256 * 1e-9;
		%20 AGC_1 (not corrected) dB/100 4 sl (see note 3)
        burst.ATT1_science = -1.0*fread(fid,1,'int32')/100;
		%21 AGC_2 (not corrected) dB/100 4 sl (see note 3)
        burst.ATT2_science = -1.0*fread(fid,1,'int32')/100;
        burst.att_sar_ku_isp = 62-(burst.ATT1_science);
            
%         burst.h0_comp_sar_isp
%         burst.cor2_comp_sar_isp
%         burst.cai_sar_isp
%         burst.fai_sar_isp
        
		%22 Total Fixed Gain Rx 1 dB/100 4 sl (see note 3)
        burst.tot_fixed_gain_1 = fread(fid,1,'int32')/100;
		%23 Total Fixed Gain Rx 2 dB/100 4 sl (see note 3)
        burst.tot_fixed_gain_2 = fread(fid,1,'int32')/100;
		%24 Transmit Power Micro-Watts 4 sl
        burst.transmit_power = fread(fid,1,'int32')* 1e-6;
        power_tx_ant_ku_chd  = mean(burst.transmit_power);
		%25 Doppler range correction (Radial component)mm 4 sl
        burst.doppler_range_correction = fread(fid,1,'int32') * 1e-3;
		%26 Instrument Range Correction tx-rx antenna mm 4 sl
        burst.instrument_range_correction_tx_rx = fread(fid,1,'int32') * 1e-3;
		%27 Instrument Range Correction mm 4 sl
        burst.instrument_range_correction_rx = fread(fid,1,'int32') * 1e-3;
		%28 Instrument Sigma 0 correction, tx-rx antenna dB/100 4 sl (see note
        burst.instrument_sigma0_correction_tx_rx = fread(fid,1,'int32')/100;
		%29 Instrument Sigma 0 correction rx only antenna dB/100 4 sl (see note
        burst.instrument_sigma0_correction_rx = fread(fid,1,'int32')/100;
		%30 Internal Phase Correction Microradians 4 sl
        burst.internal_phase_correction = fread(fid,1,'int32') * 1e-6;
		%31 External Phase Correction Microradians 4 sl
        burst.external_phase_correction = fread(fid,1,'int32')  * 1e-6;
		%32 Noise power measurement dB/100 4 sl (see note
        burst.noise_power = fread(fid,1,'int32') /100;
		%33 Phase Slope Correction Microradians 4 sl (see note
        burst.phase_slope_correction = fread(fid,1,'int32') * 1e-6;

        %EM 13.04.2016
        %compute the elapsed time since switch of SIRAL to include power
        %drift
%         time_sec_Ku_ref_time_SIRAL=burst.days + burst.seconds + round(burst.microseconds); % in seconds
%         absolute_date_burst=datevec(datestr(addtodate(datenum(date_ref_time_SIRAL),time_sec_Ku_ref_time_SIRAL,'second'),'mmmm dd, yyyy HH:MM:SS'),'mmmm dd, yyyy HH:MM:SS');
%         elapse_time_ref_switch_SIRAL=etime(absolute_date_burst,date_ref_switch_SIRAL); % in seconds
%        elapse_time_ref_switch_SIRAL=burst.days + burst.seconds + round(burst.microseconds)-date_ref_switch_SIRAL; % in seconds
%         burst.gain_corr_instr_sar = burst.tot_fixed_gain_1-(burst.ATT1_science+burst.ATT2_science)+...
%                                     burst.instrument_sigma0_correction_tx_rx+...
%                                     onboard_proc_sar_chd+elapse_time_ref_switch_SIRAL*PTR_power_drift_slope_sec...
%                                     -10*log10(zp_fact_range_cnf)+gain_scale_hamm;%
        %burst.gain_corr_instr_sar = burst.tot_fixed_gain_1-(burst.ATT1_science+burst.ATT2_science)+burst.instrument_sigma0_correction_tx_rx+onboard_proc_sar_chd-10*log10(zp_fact_range_cnf)+gain_scale_hamm+elapse_time_ref_switch_SIRAL*PTR_power_drift_slope_sec;%   
        gain_ground_processing=gain_scale_hamm+gain_scale_win_rg+10*log10(N_samples)+10*log10(N_ku_pulses_burst_chd)-10*log10(zp_fact_range_cnf);
        burst.gain_corr_instr_sar_1 = burst.tot_fixed_gain_1-(burst.ATT1_science+burst.ATT2_science)+burst.instrument_sigma0_correction_tx_rx+onboard_proc_sar_chd+gain_ground_processing;%   
        if strcmp(mode,'SIN')
            burst.gain_corr_instr_sar_2 = burst.tot_fixed_gain_2-(burst.ATT1_science+burst.ATT2_science)+burst.instrument_sigma0_correction_rx+onboard_proc_sar_chd+gain_ground_processing;%   
        end
        %34 spares 4*1 uc
        fread(fid,4,'uchar');

        
    end