function [lat,lon,alt,range_delay]= readnetCDF_Sentinel3_lat_lon(filename_L1A)
    global  semi_major_axis_cst flat_coeff_cst  N_bursts_cycle_sar_chd N_ku_pulses_burst_chd c_cst
    global bw_ku_chd
    x_pos = ncread(filename_L1A,'x_pos_l1a_echo_sar_ku');
    y_pos = ncread(filename_L1A,'y_pos_l1a_echo_sar_ku');
    z_pos = ncread(filename_L1A,'z_pos_l1a_echo_sar_ku');
    lla = ecef2lla([double(x_pos),double(y_pos),double(z_pos)],flat_coeff_cst,semi_major_axis_cst);
    lat = lla(:,1);
    lon = lla(:,2);
    alt                 = ncread(filename_L1A,'alt_l1a_echo_sar_ku');
%     alt_offset          = ncreadatt(filename_L1A,'alt_l1a_echo_sar_ku','add_offset');
%     alt_scaling_factor  = ncreadatt(filename_L1A,'alt_l1a_echo_sar_ku','scale_factor');
%     alt                 = double(alt).*double(alt_scaling_factor)+double(alt_offset);

    range_delay                 = ncread(filename_L1A,'range_ku_l1a_echo_sar_ku');
    

%     range_delay_offset          = ncreadatt(filename_L1A,'range_ku_l1a_echo_sar_ku','add_offset');
%     range_delay_scaling_factor  = ncreadatt(filename_L1A,'range_ku_l1a_echo_sar_ku','scale_factor');
%     range_delay                 = double(range_delay).*double(range_delay_scaling_factor)+double(range_delay_offset);

% 
%     cor2_applied_l1a_echo_sar_ku                 = ncread(filename_L1A,'cor2_applied_l1a_echo_sar_ku');
%     h0_applied_l1a_echo_sar_ku                   = ncread(filename_L1A,'h0_applied_l1a_echo_sar_ku');
%     burst_count_cycle_l1a_echo_sar_ku            = ncread(filename_L1A,'burst_count_cycle_l1a_echo_sar_ku');
%     pitch                                        = ncread(filename_L1A,'pitch_sat_pointing_l1a_echo_sar_ku');
%     
%     


%% COMPUTE window delay from H0 and COR2
% HPR_b=floor(double(cor2_applied_l1a_echo_sar_ku)/N_bursts_cycle_sar_chd);
% h = zeros(length(HPR_b),N_ku_pulses_burst_chd);
% for i_pulse=1:N_ku_pulses_burst_chd
%     h(:,i_pulse) = double(h0_applied_l1a_echo_sar_ku) + (double(burst_count_cycle_l1a_echo_sar_ku)-1) .* floor(HPR_b/2^4);
% end
% for i_record=1:length(HPR_b)
% CAI_tmp = floor(h(i_record,1) / 2^8); %[truncate result, algorithm b)]
% FAI_tmp = h(i_record,1) - CAI_tmp* 2^8;
% if (0 <= FAI_tmp)&& (FAI_tmp<= 127)
%     FAI = FAI_tmp;
%     CAI = CAI_tmp;
%     rx_delay(i_record) =  CAI *2^8*3.125/64*10^-9;
% 
% elseif (128 <= FAI_tmp)&& (FAI_tmp<= 255)
%     FAI = FAI_tmp - 256;
%     CAI = CAI_tmp + 1;
%     rx_delay(i_record) =  CAI *2^8*3.125/64*10^-9;
% else
%     CAI = 0;
%     disp(['ERROR in CAI computation for record ' num2str(burst_count_cycle_l1a_echo_sar_ku) ' FAI= ' num2str(FAI_tmp)]);
%     rx_delay(i_record) = double(h0_applied_l1a_echo_sar_ku(i_record))* 3.125/64*10^-9 + double(cor2_applied_l1a_echo_sar_ku(i_record))*3.125/1024*10^-9 * (double(burst_count_cycle_l1a_echo_sar_ku)-1)/N_bursts_cycle_sar_chd;  %[s] scaling factor 3.125/64*10^-9 from attributes units
% end
% end
% instr_delay = ncread(filename_L1A,'int_path_cor_ku_l1a_echo_sar_ku') ...
% + ncread(filename_L1A,'uso_cor_l1a_echo_sar_ku') ...
% + ncread(filename_L1A,'cog_cor_l1a_echo_sar_ku').* cos(pitch);
% % double(netCDF_L1A.data.int_path_cor_ku_l1a_echo_sar_ku)   .* netCDF_L1A.attributes.int_path_cor_ku_l1a_echo_sar_ku.scale_factor   + ...   %[m]
% %               double(netCDF_L1A.data.uso_cor_l1a_echo_sar_ku)           .* netCDF_L1A.attributes.uso_cor_l1a_echo_sar_ku.scale_factor           +...    %[m]
% %               double(netCDF_L1A.data.cog_cor_l1a_echo_sar_ku)           .* netCDF_L1A.attributes.cog_cor_l1a_echo_sar_ku.scale_factor .* (cos(pitch));                  %[m]
%           
% %           z_cog_ant_corr(i_burst,i_pulse) = x_cog_ant * sin(pitch_pre_dat(i_burst)) + z_cog_ant * cos(pitch_pre_dat(i_burst));      
% %           delta_pitch(i_burst,i_pulse) = - 2 * z_cog_ant_corr(i_burst,i_pulse) / c_cst;
%           
% range_delay_2 = rx_delay.' * c_cst/2 + instr_delay; %[m]
end