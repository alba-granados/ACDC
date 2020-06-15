%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the PRELIMINARY WINDOW DELAY as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: Compute the preliminary window delay.
% 
% ----------------------------------------------------------
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/09/2013)
%
% Version  record
% 1.0 2015/01/01    Imported code from S6


function [win_delay_sar_ku, win_delay_sar_ku_ref, win_delay_sar_ku_av, int_delay_cor_cal1,...
          T0_sar_pre_dat, pitch_pre_dat, yaw_pre_dat, cai_namb_sar,...
          fai_sar, mean_cai_fai_sar, z_cog_ant_corr]                       = pre_win_delay (h0_sar_isp, cor2_sar_isp, burst_sar_isp, ...
                                                                                        time_sar_ku_pre_dat,datation_tai, roll, pitch, yaw)

%% 0.Initialisation and Memory allocation
global x_cog_ant z_cog_ant N_total_bursts_sar_ku_isp i_ISP
global N_bursts_cycle_chd N_pulses_rc ext_delay_ground_chd  
global cai_cor2_unit_conv_chd h0_cor2_unit_conv_chd T0_h0_unit_conv_chd
global c_cst N_pulses_burst N_total_bursts_sar_ku
global N_total_rc int_delay_ground_chd 

N_total_rc = N_total_bursts_sar_ku_isp(i_ISP)/N_bursts_cycle_chd;

h_sar               = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
cai_namb_sar        = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
fai_sar             = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
rx_delay            = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
z_cog_ant_corr      = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
delta_pitch         = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
rx_delay_sar        = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
win_delay_sar_ku    = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
rx_delay_sar_ku_mean    = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));
instr_delay             = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));
win_delay_sar_ku_av     = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));
win_delay_sar_ku_ref    = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));
mean_cai_fai_sar        = zeros(1,N_total_bursts_sar_ku_isp(i_ISP));

cor2n           = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);
cor2_inc        = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_pulses_burst);

%% 0. Selections
[~,pitch_pre_dat,yaw_pre_dat]           = attitude_selection(time_sar_ku_pre_dat,datation_tai, roll, pitch, yaw);
T0_sar_pre_dat                          = uso_selection(time_sar_ku_pre_dat); %Currently, it returns a vector of T0_nom!
int_delay_cor_cal1                      = cal1_selection(time_sar_ku_pre_dat); %Currently, it returns a vector of 0s!

int_delay_cor_cal1 = int_delay_cor_cal1 + int_delay_ground_chd;

%%
for i_burst = 1:N_total_bursts_sar_ku_isp(i_ISP)
    progressbar([],[],(i_burst)/N_total_bursts_sar_ku_isp(i_ISP),[],[],[],[],[],[],[],[],[],[],[]);
    
    for i_pulse = 1:N_pulses_burst
    %% 1.Computation of H0
        cor2n(i_burst,i_pulse) = ((burst_sar_isp(i_burst)-1)*N_pulses_burst + (i_pulse-1)) * cor2_sar_isp(i_burst);
        if cor2n(i_burst,i_pulse) > 0
            cor2_inc(i_burst,i_pulse) = floor(cor2n(i_burst,i_pulse)/N_pulses_rc);
        else
            cor2_inc(i_burst,i_pulse) = ceil(cor2n(i_burst,i_pulse)/N_pulses_rc);
        end
        h_sar(i_burst,i_pulse) = h0_sar_isp(i_burst) * h0_cor2_unit_conv_chd + cor2_inc(i_burst,i_pulse);
        
    %% 2.Separation of CAI and FAI
        cai_namb_sar(i_burst,i_pulse) = floor(h_sar(i_burst,i_pulse) / cai_cor2_unit_conv_chd);
        fai_sar(i_burst,i_pulse) = h_sar(i_burst,i_pulse) - cai_namb_sar(i_burst,i_pulse)* cai_cor2_unit_conv_chd;
        
    %% 3.Computation of Rx delay
        rx_delay(i_burst,i_pulse) = h_sar(i_burst,i_pulse)*T0_sar_pre_dat(i_burst)/T0_h0_unit_conv_chd/h0_cor2_unit_conv_chd;
        
    %% 4. Attitude corrections
        z_cog_ant_corr(i_burst,i_pulse) = x_cog_ant * sin(pitch_pre_dat(i_burst)) + z_cog_ant * cos(pitch_pre_dat(i_burst));
        delta_pitch(i_burst,i_pulse) = - 2 * z_cog_ant_corr(i_burst,i_pulse) / c_cst;
        rx_delay_sar(i_burst,i_pulse) = rx_delay(i_burst,i_pulse) + delta_pitch(i_burst,i_pulse);
    end
    
    %% 5.Averaging
    rx_delay_sar_ku_mean(i_burst) = mean(rx_delay_sar(i_burst,2:N_pulses_burst));
    
    %% 6.Instrument corrections
    instr_delay(i_burst) = int_delay_cor_cal1(i_burst) + ext_delay_ground_chd;

    %% 7.Final window delay
    win_delay_sar_ku_av(i_burst) = rx_delay_sar_ku_mean(i_burst) + instr_delay(i_burst);
    win_delay_sar_ku_ref(i_burst) = rx_delay_sar(i_burst,2) + instr_delay(i_burst);
    for i_pulse = 2:N_pulses_burst
        win_delay_sar_ku(i_burst,i_pulse-1) = rx_delay_sar(i_burst,i_pulse) + instr_delay(i_burst);
    end
    



% % % % %% Mean CAI/FAI ==> REVIEW HOW IT IS COMPUTED!! better to sum the
% pulses during the 'for' and here only perform the mean of the 7 bursts in
% a radary cycle

    mean_cai_fai_sar(i_burst) = mean(cai_namb_sar(i_burst)*cai_cor2_unit_conv_chd + fai_sar(i_burst));
end

end
