%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% GPP 
% This code creates new burst locations and time for a given init point and
% the direction
%
% 
% ----------------------------------------------------------
% Author:    Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/09/2013)
%
% Version  record
% 1.0 2016/03/17    Imported code from S6

function [L1A] = create_new_burts_locations(L1A, direction, number_of_bursts)
global bri_nom flat_coeff_cst semi_major_axis_cst
global N_ku_pulses_burst_chd N_samples_sar_chd

if(strcmp(direction,'back'))
    time_sar_ku_gap = L1A.time_sar_ku(1)    - bri_nom(1) .* (number_of_bursts:-1:1);
    win_delay_sar_ku_gap = L1A.win_delay_sar_ku(1) * ones(1, number_of_bursts);
    reference_burst = 1;
    start_burst = reference_burst;
    end_burst = reference_burst + number_of_bursts;
elseif(strcmp(direction,'fore'))
	time_sar_ku_gap = L1A.time_sar_ku(end)  + bri_nom(1) * (1:number_of_bursts);
    win_delay_sar_ku_gap = L1A.win_delay_sar_ku(end) * ones(1, number_of_bursts);
    reference_burst = length(L1A.win_delay_sar_ku);
    start_burst = reference_burst - number_of_bursts;
    end_burst = reference_burst;
end
    
    % Read from GAP positions OSV

        % OSV file not available
        x_sar_sat_gap = interp1 (L1A.time_sar_ku(start_burst:end_burst), L1A.x_sar_sat(start_burst:end_burst), time_sar_ku_gap, 'linear', 'extrap');
        y_sar_sat_gap = interp1 (L1A.time_sar_ku(start_burst:end_burst), L1A.y_sar_sat(start_burst:end_burst), time_sar_ku_gap, 'linear', 'extrap');
        z_sar_sat_gap = interp1 (L1A.time_sar_ku(start_burst:end_burst), L1A.z_sar_sat(start_burst:end_burst), time_sar_ku_gap, 'linear', 'extrap');
        alt_rate_sar_sat_gap = interp1 (L1A.time_sar_ku(start_burst:end_burst),L1A.alt_rate_sar_sat(start_burst:end_burst), time_sar_ku_gap, 'linear', 'extrap');
        x_vel_sat_sar_gap = interp1 (L1A.time_sar_ku(start_burst:end_burst), L1A.x_vel_sat_sar(start_burst:end_burst), time_sar_ku_gap, 'linear', 'extrap');
        y_vel_sat_sar_gap = interp1 (L1A.time_sar_ku(start_burst:end_burst), L1A.y_vel_sat_sar(start_burst:end_burst), time_sar_ku_gap, 'linear', 'extrap');
        z_vel_sat_sar_gap = interp1 (L1A.time_sar_ku(start_burst:end_burst), L1A.z_vel_sat_sar(start_burst:end_burst), time_sar_ku_gap, 'linear', 'extrap');
        
        
        lla = ecef2lla([x_sar_sat_gap',y_sar_sat_gap',z_sar_sat_gap'],flat_coeff_cst,semi_major_axis_cst);
        lat_sat_gap = lla(:,1).';
        lon_sat_gap = lla(:,2).';
        alt_sat_gap = lla(:,3).';
       

        % OSV file available
%     [x_sar_sat_gap, y_sar_sat_gap, z_sar_sat_gap,...
%      x_vel_sat_sar_gap, y_vel_sat_sar_gap, z_vel_sat_sar_gap,...
%      lat_sat_gap, lon_sat_gap, alt_sat_gap,...
%      alt_rate_sar_sat_gap,~]...
%         = osv_selection(number_of_bursts,time_sar_ku_gap,1,1,1);
    
    % Interpolate GAP attitude from the att selected
    roll_gap    = zeros(1,number_of_bursts)+ L1A.roll_sar(reference_burst);
    pitch_gap   = zeros(1,number_of_bursts)+ L1A.pitch_sar(reference_burst);
    yaw_gap     = zeros(1,number_of_bursts)+ L1A.yaw_sar(reference_burst);
    pri_sar_gap = zeros(1,number_of_bursts)+ L1A.pri_sar(reference_burst);
    T0_sar_gap  = zeros(1,number_of_bursts)+ L1A.T0_sar(reference_burst);
    
    %the processID for the gap bursts is set to 0 in order to flag them and update the stack mask
    Process_ID_gap = zeros(1,number_of_bursts);
    
        % Concatenate GAP between ISPs by initialisating them
    L1A.win_delay_sar_ku    = concatenate_where (L1A.win_delay_sar_ku   , win_delay_sar_ku_gap  ,direction);
    L1A.time_sar_ku         = concatenate_where (L1A.time_sar_ku        , time_sar_ku_gap       ,direction);
    L1A.x_sar_sat           = concatenate_where (L1A.x_sar_sat          , x_sar_sat_gap         ,direction);
    L1A.y_sar_sat           = concatenate_where (L1A.y_sar_sat          , y_sar_sat_gap         ,direction);
    L1A.z_sar_sat           = concatenate_where (L1A.z_sar_sat          , z_sar_sat_gap         ,direction);

    L1A.alt_rate_sar_sat    = concatenate_where (L1A.alt_rate_sar_sat   , alt_rate_sar_sat_gap  ,direction);
    L1A.lat_sar_sat         = concatenate_where (L1A.lat_sar_sat        , lat_sat_gap           ,direction);
    L1A.lon_sar_sat         = concatenate_where (L1A.lon_sar_sat        , lon_sat_gap           ,direction);
    L1A.alt_sar_sat         = concatenate_where (L1A.alt_sar_sat        , alt_sat_gap           ,direction);
    
    L1A.x_vel_sat_sar       = concatenate_where (L1A.x_vel_sat_sar      , x_vel_sat_sar_gap     ,direction);
    L1A.y_vel_sat_sar       = concatenate_where (L1A.y_vel_sat_sar      , y_vel_sat_sar_gap     ,direction);
    L1A.z_vel_sat_sar       = concatenate_where (L1A.z_vel_sat_sar      , z_vel_sat_sar_gap     ,direction);
    L1A.pitch_sar           = concatenate_where (L1A.pitch_sar          , pitch_gap             ,direction);
    L1A.roll_sar            = concatenate_where (L1A.roll_sar           , roll_gap              ,direction);
    L1A.yaw_sar             = concatenate_where (L1A.yaw_sar            , yaw_gap               ,direction);
    L1A.pri_sar             = concatenate_where (L1A.pri_sar            , pri_sar_gap           ,direction);
    L1A.T0_sar              = concatenate_where (L1A.T0_sar             , T0_sar_gap            ,direction);
    L1A.ProcessID          = concatenate_where (L1A.ProcessID         , Process_ID_gap        ,direction);

    L1A.wfm_cal_gain_corrected = concatenate_where_3D(L1A.wfm_cal_gain_corrected, zeros(number_of_bursts,N_ku_pulses_burst_chd,N_samples_sar_chd) ,direction);
    
    L1A.N_total_bursts_sar_ku = L1A.N_total_bursts_sar_ku + number_of_bursts;

     

end
 

