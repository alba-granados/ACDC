% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT Ltd. 
% ---------------------------------------------------------
% DeDop 
% ---------------------------------------------------------
% Objective: Check if there ara data gaps on the file
% 
% ----------------------------------------------------------
% Author: Albert Garcia / isardSAT
%
% Version  record
% v1.0 first version
% v1.1 bri_nom changed for bri_sar_nom

function [L1A,original_burst_index, N_bursts] = check_continuity (last_valid_L1A,L1A,original_burst_index,i_burst,N_bursts)

global bri_sar_nom flat_coeff_cst semi_major_axis_cst c_cst

%check current time with last time


if((L1A.time_sar_ku-last_valid_L1A.time_sar_ku)> bri_sar_nom(1)*1.5)
    
    %copy L1A as the next_valid_burst to use it in the interpolations
    new_valid_L1A = L1A;
    
    %reset the current burst info to create the missing burst.
    [L1A] = create_L1A_struct;

    %fill it with new values
    L1A.time_sar_ku = last_valid_L1A.time_sar_ku  + bri_sar_nom(1);
    L1A.win_delay_sar_ku = interp1 ([last_valid_L1A.time_sar_ku, new_valid_L1A.time_sar_ku], [last_valid_L1A.win_delay_sar_ku, new_valid_L1A.win_delay_sar_ku], L1A.time_sar_ku, 'linear'); 
    L1A.x_sar_sat = interp1 ([last_valid_L1A.time_sar_ku, new_valid_L1A.time_sar_ku], [last_valid_L1A.x_sar_sat, new_valid_L1A.x_sar_sat], L1A.time_sar_ku, 'linear');
    L1A.y_sar_sat = interp1 ([last_valid_L1A.time_sar_ku, new_valid_L1A.time_sar_ku], [last_valid_L1A.y_sar_sat, new_valid_L1A.y_sar_sat], L1A.time_sar_ku, 'linear');
    L1A.z_sar_sat = interp1 ([last_valid_L1A.time_sar_ku, new_valid_L1A.time_sar_ku], [last_valid_L1A.z_sar_sat, new_valid_L1A.z_sar_sat], L1A.time_sar_ku, 'linear');
    L1A.alt_rate_sar_sat = interp1 ([last_valid_L1A.time_sar_ku, new_valid_L1A.time_sar_ku],[last_valid_L1A.alt_rate_sar_sat, new_valid_L1A.alt_rate_sar_sat], L1A.time_sar_ku, 'linear');
    L1A.x_vel_sat_sar = interp1 ([last_valid_L1A.time_sar_ku, new_valid_L1A.time_sar_ku], [last_valid_L1A.x_vel_sat_sar, new_valid_L1A.x_vel_sat_sar], L1A.time_sar_ku, 'linear');
    L1A.y_vel_sat_sar = interp1 ([last_valid_L1A.time_sar_ku, new_valid_L1A.time_sar_ku], [last_valid_L1A.y_vel_sat_sar, new_valid_L1A.y_vel_sat_sar], L1A.time_sar_ku, 'linear');
    L1A.z_vel_sat_sar = interp1 ([last_valid_L1A.time_sar_ku, new_valid_L1A.time_sar_ku], [last_valid_L1A.z_vel_sat_sar, new_valid_L1A.z_vel_sat_sar], L1A.time_sar_ku, 'linear');

    lla = ecef2lla([L1A.x_sar_sat',L1A.y_sar_sat',L1A.z_sar_sat'],flat_coeff_cst,semi_major_axis_cst);
    L1A.lat_sar_sat = lla(:,1).';
    L1A.lon_sar_sat = lla(:,2).';
    L1A.alt_sar_sat = lla(:,3).';

    L1A.roll_sar    =  last_valid_L1A.roll_sar;
    L1A.pitch_sar   =  last_valid_L1A.pitch_sar;
    L1A.yaw_sar     =  last_valid_L1A.yaw_sar;
    L1A.pri_sar     =  last_valid_L1A.pri_sar;
    L1A.T0_sar      =  last_valid_L1A.T0_sar;

    L1A.lat_sar_surf = L1A.lat_sar_sat;
    L1A.lon_sar_surf = L1A.lon_sar_sat;
    L1A.alt_sar_surf = L1A.alt_sar_sat - L1A.win_delay_sar_ku * c_cst/2;

    % geod2cart(SURF)
    p = lla2ecef([L1A.lat_sar_surf,L1A.lon_sar_surf,L1A.alt_sar_surf],flat_coeff_cst,semi_major_axis_cst);
    L1A.x_sar_surf = p(1).';
    L1A.y_sar_surf = p(2).';
    L1A.z_sar_surf = p(3).';
    L1A.confi_block_degraded=2;
    [~,L1A.doppler_ang_sar_sat] = compute_height_rate(1, L1A.x_vel_sat_sar, L1A.y_vel_sat_sar, L1A.z_vel_sat_sar,L1A.x_sar_sat ,L1A.y_sar_sat ,L1A.z_sar_sat ,L1A.x_sar_surf,L1A.y_sar_surf,L1A.z_sar_surf);

    
    % update original_burst_index,i_burst,N_bursts to read the first valid
    % burst in the future
    original_burst_index= [original_burst_index(1:i_burst-1),0,original_burst_index(i_burst:N_bursts)];
    N_bursts = N_bursts +1 ;
    
    
end



end

