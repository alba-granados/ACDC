% 20170329 v1.1 read_L1A_record function calling corrected passing N_bursts_original to it 

function  [L1A,end_of_file] = fill_empty_burst(files, last_valid_L1A, L1A,original_burst_index,N_bursts_original,i_burst)

global bri_nom flat_coeff_cst semi_major_axis_cst c_cst
% global N_ku_pulses_burst_chd N_samples_sar_chd
end_of_file 	= 0;
degraded_burst=1;

while(degraded_burst)
        original_burst_index = original_burst_index+1;
        if(original_burst_index>N_bursts_original)
            end_of_file=1;
            return;
        end
        [new_valid_L1A,files] = read_L1A_record(files,L1A,original_burst_index,i_burst,N_bursts_original);
        if(new_valid_L1A.confi_block_degraded==0)
            degraded_burst = 0;  
        end
end

% interpolate the missing burst between last_valid_L1A and new_valid_L1A.

L1A.time_sar_ku = last_valid_L1A.time_sar_ku  + bri_nom(1);
L1A.win_delay_sar_ku = last_valid_L1A.win_delay_sar_ku ;
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

[~,L1A.doppler_ang_sar_sat] = compute_height_rate(1, L1A.x_vel_sat_sar, L1A.y_vel_sat_sar, L1A.z_vel_sat_sar,L1A.x_sar_sat ,L1A.y_sar_sat ,L1A.z_sar_sat ,L1A.x_sar_surf,L1A.y_sar_surf,L1A.z_sar_surf);



end
    