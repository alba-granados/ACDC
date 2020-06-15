% Adapt S3 to S6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% S3 MPC 
% This code implements the Adaptation of the S3 L0 into the S6 L0
% INPUTs: 
%
%     Name                      Units	Type	Origin
%     S3 netCDF structure         
% 
% OUTPUTs:
%     Name                      Units	Type	Destination
%     L0 structure


function [data]     = adapt_S3L0(netCDF)


data.source_seq_count_sar_isp   = netCDF.data.seq_count_lo_echo_sar_ku;         %uint16     / -
data.days                       = netCDF.data.UTC_day_lo_echo_sar_ku;           %int16      / day since 100-01-01
data.seconds                    = netCDF.data.UTC_sec_lo_echo_sar_ku;           %float64    / seconds in the day
data.microseconds               = netCDF.data.sral_fine_time_lo_echo_sar_ku;    %uint32     / 137.5 5*10e-9 seconds
modeID_mode
data.ProcessID

switch modeID_mode
    case 1 %LRM
        disp('LRM product');
        return
    case 2 %SAR
       data.ProcessID(i_burst) = 58;

end

data.inst_id_sar_isp            = netCDF.data.oper_instr_lo_echo_sar_ku;        %int8       / -
% data.pri_sar_isp(i_burst) = pri_sar_chd;
% data.ambiguity_order_sar_isp  = 0;
data.burst_sar_ku               = netCDF.data.burst_nb_lo_echo_sar_ku;          %int8       / count
data.burst_sar_ku_fbr           = 1+mod(data.burst_sar_ku,N_bursts_cycle_chd);
data.lat_sar_sat                = netCDF.data.lat_lo_echo_sar_ku;               %int32      / degreer north               
data.lon_sar_sat                = netCDF.data.lon_lo_echo_sar_ku;               %int32      / degreer east
data.alt_sar_sat                = netCDF.data.alt_lo_echo_sar_ku;               %int32      / m
data.alt_rate_sar_sat           = netCDF.data.alt_rate_lo_echo_sar_ku;          %int16      / m/s
data.x_vel_sat_sar              = netCDF.data.x_vel_lo_echo_sar_ku;             %float64    / m/s
data.y_vel_sat_sar              = netCDF.data.y_vel_lo_echo_sar_ku;             %float64    / m/s
data.z_vel_sat_sar              = netCDF.data.z_vel_lo_echo_sar_ku;             %float64    / m/s

% data.roll_sar % [rad]
% data.pitch_sar % [rad]
% data.yaw_sar % [rad]


% %15 Window Delay (2way) uncorrected for instrument delays 10-12 s 8 sll
% data.win_delay_sar_ku

data.h0_comp_sar_isp            = netCDF.data.ho_nav_dem_lo_echo_sar_ku;        %uint32     / 3.125/64*10e-9 s
data.h0_sar_isp                 = netCDF.data.ho_applied_lo_echo_sar_ku;        %uint32     / 3.125/64*10e-9 s
data.cor2_comp_sar_isp          = netCDF.data.cor2_nav_dem_lo_echo_sar_ku;      %uint32     / 3.125/1024*10e-9 s
data.cor2_sar_isp               = netCDF.data.cor2_applied_lo_echo_sar_ku;      %uint32     / 3.125/1024*10e-9 s
data.ATT1_science               = netCDF.data.agccode_ku_lo_echo_sar_ku;        %int8       / dB
data.ATT2_science
data.att_sar_ku_isp = 62-(data.ATT1_science);

%% CALIBRATION GROUP to be filled during L0-L1A
% data.tot_fixed_gain_1
% data.tot_fixed_gain_2
% data.transmit_power
% data.doppler_range_correction
% data.instrument_range_correction_tx_rx
% data.instrument_range_correction
% data.instrument_sigma0_correction_tx_rx
% 
% data.internal_phase_correction
% data.external_phase_correction
% data.noise_power
% data.phase_slope_correction

% data.nimp_sar_isp

wfm_iq_sar_ku_i         = netCDF.data.i_meas_lo_echo_sar_ku;                %int8       / -
wfm_iq_sar_ku_q         = netCDF.data.q_meas_lo_echo_sar_ku;                %int8       / -

data.wfm_iq_sar_ku_isp     = (wfm_iq_sar_ku_i + 1i * wfm_iq_sar_ku_q);





end