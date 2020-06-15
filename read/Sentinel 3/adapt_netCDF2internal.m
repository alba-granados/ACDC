%v1.0 First Version
%v1.1 Added Calibration (Power, CoG to antenna)
%v1.2 Added new cnf variables move_CoG2ant_alongtrack_cnf window_delay_source_cnf
function [L1A] = adapt_netCDF2internal(netCDF_L1A,filename_L1A, i_burst)


global N_ku_pulses_burst_chd N_bursts_cycle_sar_chd
global mode N_samples 
global flat_coeff_cst semi_major_axis_cst
global pri_sar_chd bri_sar_chd prf_sar_chd brf_sar_chd T0_chd
global pri_sar_nom bri_sar_nom prf_sar_nom brf_sar_nom c_cst
global height_rate_application_cnf window_delay_source_cnf move_CoG2ant_alongtrack_cnf
% global diff_wd_delay instr_delay_g
L1A         = [];
eccentricity = sqrt(flat_coeff_cst*(2-flat_coeff_cst));
L1A = setfield(L1A,'source_seq_count_sar_isp',netCDF_L1A.data.seq_count_l1a_echo_sar_ku);
L1A = setfield(L1A,'days',double(netCDF_L1A.data.UTC_day_l1a_echo_sar_ku)); 
L1A = setfield(L1A,'seconds',double(netCDF_L1A.data.UTC_sec_l1a_echo_sar_ku));
L1A = setfield(L1A,'microseconds',0);
L1A = setfield(L1A,'USO_correction',double(netCDF_L1A.data.uso_cor_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.uso_cor_l1a_echo_sar_ku.scale_factor));
L1A = setfield(L1A,'ProcessID',58);
L1A = setfield(L1A,'source_seq_count_sar_ku_fbr',0);
L1A = setfield(L1A,'ins_id',0);
L1A = setfield(L1A,'ins_tracking_mode',0);
L1A = setfield(L1A,'ins_loop_stat',0);
L1A = setfield(L1A,'inst_id_sar_isp',0);
L1A = setfield(L1A,'pri_sar_isp',pri_sar_chd);
L1A = setfield(L1A,'ambiguity_order_sar_isp',0);
L1A = setfield(L1A,'burst_sar_ku',netCDF_L1A.data.burst_count_prod_l1a_echo_sar_ku);
L1A = setfield(L1A,'burst_sar_ku_fbr',netCDF_L1A.data.burst_count_prod_l1a_echo_sar_ku);

pitch = double(netCDF_L1A.data.pitch_sat_pointing_l1a_echo_sar_ku).*netCDF_L1A.attributes.pitch_sat_pointing_l1a_echo_sar_ku.scale_factor;
yaw = double(netCDF_L1A.data.yaw_sat_pointing_l1a_echo_sar_ku).*netCDF_L1A.attributes.yaw_sat_pointing_l1a_echo_sar_ku.scale_factor;
x_cog_ant= 0.683; % [meters] Provided by François Boy FIXED VALUE as it is not in the L1A product
z_cog_ant=double(netCDF_L1A.data.cog_cor_l1a_echo_sar_ku)           .* netCDF_L1A.attributes.cog_cor_l1a_echo_sar_ku.scale_factor; 

if(move_CoG2ant_alongtrack_cnf)

    try
        x_pos_pair = ncread(filename_L1A,'x_pos_l1a_echo_sar_ku',i_burst,2);
        y_pos_pair = ncread(filename_L1A,'y_pos_l1a_echo_sar_ku',i_burst,2);
        z_pos_pair = ncread(filename_L1A,'z_pos_l1a_echo_sar_ku',i_burst,2);



        x_cog_ant_aux = x_cog_ant * cos(pitch) + z_cog_ant * sin(pitch);
        x_cog_ant_corr = x_cog_ant_aux * cos(yaw);
        vel=sqrt((double(netCDF_L1A.data.x_vel_l1a_echo_sar_ku))^2+(double(netCDF_L1A.data.y_vel_l1a_echo_sar_ku))^2+(double(netCDF_L1A.data.z_vel_l1a_echo_sar_ku))^2);
        % % vlla = ecef2lla([double(netCDF_L1A.data.x_vel_l1a_echo_sar_ku),double(netCDF_L1A.data.y_vel_l1a_echo_sar_ku),double(netCDF_L1A.data.z_vel_l1a_echo_sar_ku)],flat_coeff_cst,semi_major_axis_cst);
        % vlat_sat = vlla(1);
        % % vlon_sat = vlla(2);
        % % vh_sat = vlla(3);
        t_cog_ant_corr = x_cog_ant_corr / vel;
        corrected_time = double(netCDF_L1A.data.time_l1a_echo_sar_ku)+t_cog_ant_corr;
    catch
        x_cog_ant_corr= x_cog_ant;
        x_pos_pair = ncread(filename_L1A,'x_pos_l1a_echo_sar_ku',i_burst-1,2);
        y_pos_pair = ncread(filename_L1A,'y_pos_l1a_echo_sar_ku',i_burst-1,2);
        z_pos_pair = ncread(filename_L1A,'z_pos_l1a_echo_sar_ku',i_burst-1,2);
    end

    lla = ecef2lla([x_pos_pair,y_pos_pair,z_pos_pair],flat_coeff_cst,semi_major_axis_cst); % error when computing the last one

    [arclen_v,az_v] = distance('gc',lla(1,1),lla(1,2),lla(2,1),lla(2,2),[semi_major_axis_cst eccentricity]);
    [lats,lons] = reckon(lla(1,1),lla(1,2),x_cog_ant_corr,az_v,[semi_major_axis_cst eccentricity]);
    lat_sat=lats;
    lon_sat=lons;
else
    
    lla = ecef2lla([double(netCDF_L1A.data.x_pos_l1a_echo_sar_ku),double(netCDF_L1A.data.y_pos_l1a_echo_sar_ku),double(netCDF_L1A.data.z_pos_l1a_echo_sar_ku)],flat_coeff_cst,semi_major_axis_cst);
    lat_sat = lla(1,1);
    lon_sat = lla(1,2);
    
end


% lat_sat = lla(1,1);
% lon_sat = lla(1,2);

alt_sat2 = lla(1,3); 
L1A = setfield(L1A,'lat_sar_sat',lat_sat);
L1A = setfield(L1A,'lon_sar_sat',lon_sat);
alt_sat = double(netCDF_L1A.data.alt_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.alt_l1a_echo_sar_ku.scale_factor)+netCDF_L1A.attributes.alt_l1a_echo_sar_ku.add_offset;
L1A = setfield(L1A,'alt_sar_sat',alt_sat);
%the difference using the alt_sat instead of alt_sat2 is -4.5009e-05 meters for the record cheched

alt_rate=double(netCDF_L1A.data.orb_alt_rate_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.orb_alt_rate_l1a_echo_sar_ku.scale_factor)+netCDF_L1A.attributes.orb_alt_rate_l1a_echo_sar_ku.add_offset;
L1A = setfield(L1A,'alt_rate_sar_sat',alt_rate);

L1A = setfield(L1A,'x_vel_sat_sar',double(netCDF_L1A.data.x_vel_l1a_echo_sar_ku));
L1A = setfield(L1A,'y_vel_sat_sar',double(netCDF_L1A.data.y_vel_l1a_echo_sar_ku));
L1A = setfield(L1A,'z_vel_sat_sar',double(netCDF_L1A.data.z_vel_l1a_echo_sar_ku));
L1A = setfield(L1A,'roll_sar',double(netCDF_L1A.data.roll_sral_mispointing_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.roll_sral_mispointing_l1a_echo_sar_ku.scale_factor));
L1A = setfield(L1A,'pitch_sar',double(netCDF_L1A.data.pitch_sral_mispointing_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.pitch_sral_mispointing_l1a_echo_sar_ku.scale_factor));
L1A = setfield(L1A,'yaw_sar',double(netCDF_L1A.data.yaw_sral_mispointing_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.yaw_sral_mispointing_l1a_echo_sar_ku.scale_factor));
L1A = setfield(L1A,'confi_block_degraded',0); %double(netCDF_L1A.data.nav_bul_status_l1a_echo_sar_ku)
L1A = setfield(L1A,'mea_conf_data_sar_ku_fbr',0);


%'Distance between the altimeter reference point and the surface height 
% associated to a range gate used as reference inside the tracking window 
%(reference tracking point), corrected for USO frequency drift and internal path correction
if window_delay_source_cnf == 0 %L0
    HPR_b=floor(double(netCDF_L1A.data.cor2_applied_l1a_echo_sar_ku)/N_bursts_cycle_sar_chd);
    h = zeros(1,N_ku_pulses_burst_chd);
    for i_pulse=1:N_ku_pulses_burst_chd
        h(i_pulse) = double(netCDF_L1A.data.h0_applied_l1a_echo_sar_ku) + (double(netCDF_L1A.data.burst_count_cycle_l1a_echo_sar_ku)-1) * floor(HPR_b/2^4);
    end

    CAI_tmp = floor(h(1) / 2^8); %[truncate result, algorithm b)]
    FAI_tmp = h(1) - CAI_tmp* 2^8;
    if (0 <= FAI_tmp)&& (FAI_tmp<= 127)
        FAI = FAI_tmp;
        CAI = CAI_tmp;
        rx_delay =  CAI *2^8*3.125/64*10^-9;

    elseif (128 <= FAI_tmp)&& (FAI_tmp<= 255)
        FAI = FAI_tmp - 256;
        CAI = CAI_tmp + 1;
        rx_delay =  CAI *2^8*3.125/64*10^-9;
    else
        CAI = 0;
        disp(['ERROR in CAI computation for record ' num2str(netCDF_L1A.data.burst_count_cycle_l1a_echo_sar_ku) ' FAI= ' num2str(FAI_tmp)]);
        rx_delay = double(netCDF_L1A.data.h0_applied_l1a_echo_sar_ku)* 3.125/64*10^-9 + double(netCDF_L1A.data.cor2_applied_l1a_echo_sar_ku)*3.125/1024*10^-9 * (double(netCDF_L1A.data.burst_count_cycle_l1a_echo_sar_ku)-1)/N_bursts_cycle_sar_chd;  %[s] scaling factor 3.125/64*10^-9 from attributes units
    end
    instr_delay = double(netCDF_L1A.data.int_path_cor_ku_l1a_echo_sar_ku)   .* netCDF_L1A.attributes.int_path_cor_ku_l1a_echo_sar_ku.scale_factor   + ...   %[m]
                  double(netCDF_L1A.data.uso_cor_l1a_echo_sar_ku)           .* netCDF_L1A.attributes.uso_cor_l1a_echo_sar_ku.scale_factor +...    %[m]
                  double(netCDF_L1A.data.cog_cor_l1a_echo_sar_ku)           .* netCDF_L1A.attributes.cog_cor_l1a_echo_sar_ku.scale_factor .* (cos(pitch));                  %[m]

    %           z_cog_ant_corr(i_burst,i_pulse) = x_cog_ant * sin(pitch_pre_dat(i_burst)) + z_cog_ant * cos(pitch_pre_dat(i_burst));      
    %           delta_pitch(i_burst,i_pulse) = - 2 * z_cog_ant_corr(i_burst,i_pulse) / c_cst;

    range_delay = rx_delay * c_cst/2 + instr_delay; %[m]
    % diff_wd_delay = [diff_wd_delay range_delay_2-range_delay];
    % instr_delay_g = [instr_delay_g ; double(netCDF_L1A.data.int_path_cor_ku_l1a_echo_sar_ku)   .* netCDF_L1A.attributes.int_path_cor_ku_l1a_echo_sar_ku.scale_factor   ...   %[m]
    %               double(netCDF_L1A.data.uso_cor_l1a_echo_sar_ku)           .* netCDF_L1A.attributes.uso_cor_l1a_echo_sar_ku.scale_factor           ...    %[m]
    %               double(netCDF_L1A.data.cog_cor_l1a_echo_sar_ku)           .* netCDF_L1A.attributes.cog_cor_l1a_echo_sar_ku.scale_factor];

else %L1A
    range_delay = double(netCDF_L1A.data.range_ku_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.range_ku_l1a_echo_sar_ku.scale_factor)+netCDF_L1A.attributes.range_ku_l1a_echo_sar_ku.add_offset;
    instr_delay = double(netCDF_L1A.data.cog_cor_l1a_echo_sar_ku)           .* netCDF_L1A.attributes.cog_cor_l1a_echo_sar_ku.scale_factor .* (cos(pitch));
    range_delay = range_delay + instr_delay;
end
    L1A = setfield(L1A,'win_delay_sar_ku',range_delay/c_cst*2);

L1A = setfield(L1A,'h0_comp_sar_isp',double(netCDF_L1A.data.h0_nav_dem_l1a_echo_sar_ku)*3.125/64*10^-9);%[seconds]
L1A = setfield(L1A,'cor2_comp_sar_isp',double(netCDF_L1A.data.cor2_nav_dem_l1a_echo_sar_ku)*3.125/1024*10^-9);%[seconds]
L1A = setfield(L1A,'cai_sar_isp',double(netCDF_L1A.data.h0_applied_l1a_echo_sar_ku)*3.125/64*10^-9);%[seconds]
L1A = setfield(L1A,'fai_sar_isp',double(netCDF_L1A.data.cor2_applied_l1a_echo_sar_ku)*3.125/1024*10^-9);
L1A = setfield(L1A,'ATT1_science',-1.0*double(netCDF_L1A.data.agc_ku_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.agc_ku_l1a_echo_sar_ku.scale_factor)+netCDF_L1A.attributes.agc_ku_l1a_echo_sar_ku.add_offset);
L1A = setfield(L1A,'ATT2_science',0);
L1A = setfield(L1A,'att_sar_ku_isp',62-L1A.ATT1_science); %as done in CR2
L1A = setfield(L1A,'tot_fixed_gain_1',0);
L1A = setfield(L1A,'tot_fixed_gain_2',0);
L1A = setfield(L1A,'transmit_power',0);
L1A = setfield(L1A,'doppler_range_correction',0);
L1A = setfield(L1A,'instrument_range_correction_tx_rx',double(netCDF_L1A.data.int_path_cor_ku_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.int_path_cor_ku_l1a_echo_sar_ku.scale_factor)+netCDF_L1A.attributes.int_path_cor_ku_l1a_echo_sar_ku.add_offset);
L1A = setfield(L1A,'instrument_range_correction_rx',0);
L1A = setfield(L1A,'instrument_sigma0_correction_tx_rx',double(netCDF_L1A.data.sig0_cal_ku_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.sig0_cal_ku_l1a_echo_sar_ku.scale_factor)+netCDF_L1A.attributes.sig0_cal_ku_l1a_echo_sar_ku.add_offset); %SIGMA0 CAL
L1A = setfield(L1A,'instrument_sigma0_correction_rx',double(netCDF_L1A.data.scale_factor_ku_l1a_echo_sar_ku)*double(netCDF_L1A.attributes.scale_factor_ku_l1a_echo_sar_ku.scale_factor)+netCDF_L1A.attributes.scale_factor_ku_l1a_echo_sar_ku.add_offset); %[dB] %SCALE FACTOR SIGMA0 , GOES TO SCALING FACTOR ALGORITHM
L1A = setfield(L1A,'internal_phase_correction',0);
L1A = setfield(L1A,'external_phase_correction',0);
L1A = setfield(L1A,'noise_power',0);
L1A = setfield(L1A,'phase_slope_correction',0);
L1A = setfield(L1A,'gain_corr_instr_sar_1',0);


L1A = setfield(L1A,'dry_tropo_correction_bursts',0);
L1A = setfield(L1A,'wet_tropo_correction_bursts',0);
L1A = setfield(L1A,'inverse_baro_correction_bursts',0);
L1A = setfield(L1A,'Dynamic_atmospheric_correction_bursts',0);
L1A = setfield(L1A,'GIM_iono_correction_bursts',0);
L1A = setfield(L1A,'model_iono_correction_bursts',0);
L1A = setfield(L1A,'ocean_equilibrium_tide_bursts',0);
L1A = setfield(L1A,'long_period_tide_height_bursts',0);
L1A = setfield(L1A,'ocean_loading_tide_bursts',0);
L1A = setfield(L1A,'solid_earth_tide_bursts',0);
L1A = setfield(L1A,'geocentric_polar_tide_bursts',0);
L1A = setfield(L1A,'surface_type_flag_bursts',0);

L1A = setfield(L1A,'nimp_sar_isp',0);



i_pulse = double(netCDF_L1A.data.i_meas_ku_l1a_echo_sar_ku);
q_pulse = double(netCDF_L1A.data.q_meas_ku_l1a_echo_sar_ku);

wfm_cal_gain_UNCORRECTED = (i_pulse+1i*q_pulse).';
wfm_cal_gain_CORRECTED = waveforms_correction_S3(wfm_cal_gain_UNCORRECTED,netCDF_L1A);
L1A = setfield(L1A,'wfm_cal_gain_corrected',wfm_cal_gain_CORRECTED);

if height_rate_application_cnf
    %Assuming FAI is not applied on-board
    [L1A] = height_rate_alignment(L1A);  
end


lat_sar_surf = lat_sat;
lon_sar_surf = lon_sat;
alt_sar_surf = alt_sat - range_delay;


% geod2cart(SURF)
p = lla2ecef([lat_sar_surf,lon_sar_surf,alt_sar_surf],flat_coeff_cst,semi_major_axis_cst);
x_sar_surf = p(1).';
y_sar_surf = p(2).';
z_sar_surf = p(3).';


L1A = setfield(L1A,'time_sar_ku',double(netCDF_L1A.data.time_l1a_echo_sar_ku));
L1A = setfield(L1A,'x_sar_sat',double(netCDF_L1A.data.x_pos_l1a_echo_sar_ku));
L1A = setfield(L1A,'y_sar_sat',double(netCDF_L1A.data.y_pos_l1a_echo_sar_ku));
L1A = setfield(L1A,'z_sar_sat',double(netCDF_L1A.data.z_pos_l1a_echo_sar_ku));
L1A = setfield(L1A,'lat_sar_surf',lat_sar_surf);
L1A = setfield(L1A,'lon_sar_surf',lon_sar_surf);
L1A = setfield(L1A,'alt_sar_surf',alt_sar_surf);
L1A = setfield(L1A,'x_sar_surf',x_sar_surf);
L1A = setfield(L1A,'y_sar_surf',y_sar_surf);
L1A = setfield(L1A,'z_sar_surf',z_sar_surf);
[~,doppler_ang_sar_sat] = compute_height_rate(1, double(netCDF_L1A.data.x_vel_l1a_echo_sar_ku), double(netCDF_L1A.data.y_vel_l1a_echo_sar_ku), double(netCDF_L1A.data.z_vel_l1a_echo_sar_ku),double(netCDF_L1A.data.x_pos_l1a_echo_sar_ku) ,double(netCDF_L1A.data.y_pos_l1a_echo_sar_ku) ,double(netCDF_L1A.data.z_pos_l1a_echo_sar_ku) ,x_sar_surf,y_sar_surf,z_sar_surf);


L1A = setfield(L1A,'doppler_ang_sar_sat',doppler_ang_sar_sat);

L1A = setfield(L1A,'T0_sar',T0_chd);
L1A = setfield(L1A,'pri_sar',pri_sar_chd);

%beam angles
L1A = setfield(L1A,'N_beams_sar_ku',0);
L1A = setfield(L1A,'surf_loc_index',0);
L1A = setfield(L1A,'beam_ang',0);
L1A = setfield(L1A,'beam_ang_nadir_index',0);
L1A = setfield(L1A,'beam_ang_index',0);
L1A = setfield(L1A,'start_beam',0);
L1A = setfield(L1A,'end_beam',0);
L1A = setfield(L1A,'beam_index',0);

%azimuth processing
L1A = setfield(L1A,'beams_focused_shifted',zeros(N_ku_pulses_burst_chd,N_samples));
if(strcmp(mode,'SIN'))
    L1A = setfield(L1A,'beams_focused_shifted_2',zeros(N_ku_pulses_burst_chd,N_samples));
end
p = lla2ecef([L1A.lat_sar_sat.',L1A.lon_sar_sat.',L1A.alt_sar_sat.'],flat_coeff_cst,semi_major_axis_cst);
L1A.x_sar_sat = p(:,1).';
L1A.y_sar_sat = p(:,2).';
L1A.z_sar_sat = p(:,3).';

L1A.lat_sar_surf = L1A.lat_sar_sat;
L1A.lon_sar_surf = L1A.lon_sar_sat;
L1A.alt_sar_surf = L1A.alt_sar_sat - L1A.win_delay_sar_ku * c_cst/2;

% geod2cart(SURF)
p = lla2ecef([L1A.lat_sar_surf,L1A.lon_sar_surf,L1A.alt_sar_surf],flat_coeff_cst,semi_major_axis_cst);
L1A.x_sar_surf = p(1).';
L1A.y_sar_surf = p(2).';
L1A.z_sar_surf = p(3).';

pri_sar_nom = pri_sar_chd;
prf_sar_nom = prf_sar_chd; 
bri_sar_nom = bri_sar_chd; 
brf_sar_nom = brf_sar_chd; 


end

