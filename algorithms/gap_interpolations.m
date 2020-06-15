%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the GAP Interpolation
% algorithm 
% For the first ISP, the code just prepare the init structure to input in
% the surface locations
% ---------------------------------------------------------
% Objective: Read ISP data
% 
% INPUTs : L1A
% OUTPUTs: L1A with gaps
%
% ----------------------------------------------------------
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/02/2015)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N_total_bursts_sar_ku_isp bri_nom pri_sar_nom


% interpolate the gaps to find burst positions;
gap_bursts_time             = zeros(1,N_ISPs-1);
N_total_gap_bursts          = zeros(1,N_ISPs-1);
N_total_bursts_sar_ku_all   = zeros(1,N_ISPs-1);
i_gap=0;

for i_ISP=2:N_ISPs
    end_actual  = sum(N_total_bursts_sar_ku_isp(1:i_ISP-1))+sum(N_total_gap_bursts);
    bursts1     = init_burst(1):end_actual;
    bursts2     = end_actual+1:sum(N_total_bursts_sar_ku_isp)+sum(N_total_gap_bursts); % FIX end
    
    i_gap=i_gap+1;
    
    gap_bursts_time(i_gap) = time_sar_ku(init_burst(i_ISP)) - time_sar_ku(end_burst(i_ISP-1));
    N_total_gap_bursts(i_gap) = floor(gap_bursts_time(i_gap)/mean(bri_nom))-1;
    
    % Interpolate GAP window delays
    win_delay_sar_ku_gap = mean([win_delay_sar_ku(init_burst(i_ISP)),win_delay_sar_ku(end_burst(i_ISP-1))])*ones(1,N_total_gap_bursts(i_gap));
	
    % Find GAP times
    time_sar_ku_gap = time_sar_ku(end_burst(i_ISP-1))+mean(bri_nom)*(1:N_total_gap_bursts(i_gap));
    
    % Read from GAP positions OSV 
    [x_sar_sat_gap, y_sar_sat_gap, z_sar_sat_gap,...
     x_vel_sat_sar_gap, y_vel_sat_sar_gap, z_vel_sat_sar_gap,...
     lat_sat_gap, lon_sat_gap, alt_sat_gap,...
     alt_rate_sar_sat_gap,~]...
        = osv_selection(N_total_gap_bursts(i_gap),time_sar_ku_gap,1,1,1);
    
    % Interpolate GAP attitude from the att selected
    [roll_gap,pitch_gap,yaw_gap]...
        = attitude_selection(time_sar_ku_gap,datation_tai, roll, pitch, yaw);

    
    pri_sar_gap = zeros(1,N_total_gap_bursts(i_gap))+ mean(pri_sar_nom);
    T0_sar_gap  = zeros(1,N_total_gap_bursts(i_gap))+ mean(T0_sar);
    
    %the processID for the ga pbursts is set to 0 in order to flag them and update the stack mask
    Process_ID_gap = zeros(1,N_total_gap_bursts(i_gap)); 

    
    % Concatenate GAP between ISPs by initialisating them
    win_delay_sar_ku    = [win_delay_sar_ku(bursts1)    win_delay_sar_ku_gap    win_delay_sar_ku(bursts2)];
    time_sar_ku         = [time_sar_ku(bursts1)         time_sar_ku_gap         time_sar_ku(bursts2)];
    x_sar_sat           = [x_sar_sat(bursts1)           x_sar_sat_gap           x_sar_sat(bursts2)];
    y_sar_sat           = [y_sar_sat(bursts1)           y_sar_sat_gap           y_sar_sat(bursts2)];
    z_sar_sat           = [z_sar_sat(bursts1)           z_sar_sat_gap           z_sar_sat(bursts2)];
    alt_rate_sar_sat    = [alt_rate_sar_sat(bursts1)    alt_rate_sar_sat_gap    alt_rate_sar_sat(bursts2)];
    lat_sar_sat         = [lat_sar_sat(bursts1)         lat_sat_gap.'           lat_sar_sat(bursts2)];
    lon_sar_sat         = [lon_sar_sat(bursts1)         lon_sat_gap.'           lon_sar_sat(bursts2)];
    alt_sar_sat         = [alt_sar_sat(bursts1)         alt_sat_gap.'           alt_sar_sat(bursts2)];
  
    
    x_vel_sat_sar       = [x_vel_sat_sar(bursts1)       x_vel_sat_sar_gap   x_vel_sat_sar(bursts2)];
    y_vel_sat_sar       = [y_vel_sat_sar(bursts1)       y_vel_sat_sar_gap   y_vel_sat_sar(bursts2)];
    z_vel_sat_sar       = [z_vel_sat_sar(bursts1)       z_vel_sat_sar_gap   z_vel_sat_sar(bursts2)];
    pitch_sar           = [pitch_sar(bursts1)           pitch_gap           pitch_sar(bursts2)];
    roll_sar            = [roll_sar(bursts1)            roll_gap            roll_sar(bursts2)];
    yaw_sar             = [yaw_sar(bursts1)             yaw_gap             yaw_sar(bursts2)];
    pri_sar             = [pri_sar(bursts1)             pri_sar_gap         pri_sar(bursts2)];
    T0_sar              = [T0_sar(bursts1)              T0_sar_gap          T0_sar(bursts2)];
    %Theburst are artificially created in order to fill with zeros the beams in the stacks
    wfm_cal_gain_corrected = [  wfm_cal_gain_corrected(bursts1,:,:);...
                                zeros(N_total_gap_bursts(i_gap),N_ku_pulses_burst_chd,N_samples_sar_chd); ...
                                wfm_cal_gain_corrected(bursts2,:,:)];
                            
    Process_ID2 = [Process_ID2(bursts1), Process_ID_gap, Process_ID2(bursts2)];
    
end
N_total_bursts_sar_ku = sum(N_total_bursts_sar_ku_isp)+sum(N_total_gap_bursts);

    %>> Orbit to surface
    %     win_delay_surf = mean([actual.win_delay_surf(actual.N_total_surf_loc),next.win_delay_sar_ku(1)])*ones(1,gap.N_total_surf_loc);
    lat_surf = lat_sar_sat;
    lon_surf = lon_sar_sat;
    alt_surf = alt_sar_sat - win_delay_sar_ku * c_cst/2;

    %>> Geodetic to cartesian
    p = lla2ecef([lat_surf.',lon_surf.',alt_surf.'],flat_coeff_cst,semi_major_axis_cst);
    x_sar_surf = p(:,1).';
    y_sar_surf = p(:,2).';
    z_sar_surf = p(:,3).';
    
    
    
    
    
    
% Build the init structure to inititalisate the surface locations
init.lat_surf           = lat_surf(1);
init.lon_surf           = lon_surf(1);
init.alt_surf           = alt_surf(1);
init.x_surf             = x_sar_surf(1);
init.y_surf             = y_sar_surf(1);
init.z_surf             = z_sar_surf(1);
init.time_surf          = time_sar_ku(1);
init.lat_sat            = lat_sar_sat(1);
init.lon_sat            = lon_sar_sat(1);
init.alt_sat            = alt_sar_sat(1);
init.x_sat              = x_sar_sat (1);
init.y_sat              = y_sar_sat (1);
init.z_sat              = z_sar_sat (1);
init.x_vel_sat          = x_vel_sat_sar(1);
init.y_vel_sat          = y_vel_sat_sar(1);
init.z_vel_sat          = z_vel_sat_sar(1);
init.pitch_surf         = pitch_sar(1);
init.roll_surf          = roll_sar(1);
init.yaw_surf           = yaw_sar(1);
init.win_delay_surf     = win_delay_sar_ku(1);
init.alt_rate_sat       = alt_rate_sar_sat(1);
init.i_curr = 0;
 
% 1. Surface generation
% In order to generate the surface locations, get the coordinates from the satellite locations and the window delay:

lat_sar_surf = lat_sar_sat;
lon_sar_surf = lon_sar_sat;
alt_sar_surf = alt_sar_sat - win_delay_sar_ku * c_cst / 2;


% Then, transform the surface locations from geodetic to cartesian coordinates (see �5.3.2) and get x_surf_geoloc, y_surf_geoloc and z_surf_geoloc.
for i_burst = N_total_bursts_sar_ku
    
    p = lla2ecef([lat_sar_surf(i_burst).', lon_sar_surf(i_burst).', alt_sar_surf(i_burst).'], flat_coeff_cst, semi_major_axis_cst);
    x_sar_surf(i_burst) = p(1);
    y_sar_surf(i_burst) = p(2);
    z_sar_surf(i_burst) = p(3);
end    
% 
% 
% if(i_gap ==0)
% % If it is the first ISP, the init paramenters are the ones from the first surface    
% 
%     
% else
%     
%     % If it is the second ISP, buscarem el gap.
% 
%     
%     
%     
% %>> Surface gap
% %     gap_surf_time = data(i_ISP+1).time_surf(1) - time_surf(N_total_surf_loc(i_ISP));
% %     time_between_surfs = mean([diff(time_surf),diff(data(i_ISP+1).time_surf)]);
% %     gap_surfs = ceil(gap_surf_time/time_between_surfs);
% 
%     % NEW
%     time_between_surfs(i_gap)       = mean(diff(time_surf));
%     gap_surf_time(i_gap)            = time_sar_ku(init_burst(i_ISP)) - time_surf(end_surf(i_ISP-1));
%     gap(i_gap).N_total_surf_loc     = ceil(gap_surf_time/time_between_surfs); %We apply the 'ceil' function because we want to compute one extra surface
%     
%     
%     %handle the burst and surfaces index
%     
%     init_surf   = [init_surf init_surf(i_ISP)+gap(i_gap).N_total_surf_loc];
%     end_surf    = [end_surf end_surf(i_ISP-1)+gap(i_gap).N_total_surf_loc];
%     init_burst  = [init_burst init_burst(i_ISP)+gap(i_gap).N_total_bursts];
%     end_burst   = [end_burst(1:i_ISP-1) end_burst(i_ISP)-N_total_bursts_sar_ku_isp(i_ISP)+gap(i_gap).N_total_bursts end_burst(i_ISP)+gap(i_gap).N_total_bursts];
%     
%     
% %% Create surface locations within the gap (using a new struct, for now)
% %>> Time interpolation
% gap(i_gap).time_surf            = time_surf(end_surf(i_ISP-1)) + time_between_surfs(i_gap)*(1:gap(i_gap).N_total_surf_loc);
% 
% 
% %>> Orbit interpolation
% 
% %     N_surf_interp = 5;
% %     start_actual = N_total_surf_loc(i_ISP) - N_surf_interp;
% %     end_actual = N_total_surf_loc(i_ISP);
% %     start_next = 1;
% %     end_next = N_surf_interp;
%     
%     
%     % With this new time, call the OSV
%     [gap(i_gap).x_sat, gap(i_gap).y_sat, gap(i_gap).z_sat,...
%      gap(i_gap).x_vel_sat, gap(i_gap).y_vel_sat, gap(i_gap).z_vel_sat,...
%      gap(i_gap).lat_sat, gap(i_gap).lon_sat, gap(i_gap).alt_sat,...
%      gap(i_gap).alt_rate,~]...
%                                                 = osv_selection(gap(i_gap).N_total_surf_loc,gap(i_gap).time_surf,1,1,1);
%     
%     
% 	[gap(i_gap).roll_surf,gap(i_gap).pitch_surf,gap(i_gap).yaw_surf]...
%                                                 = attitude_selection(gap(i_gap).time_surf,datation_tai, roll, pitch, yaw);
%     
%     
% %     gap.x_surf(:) = spline( [time_surf(start_actual:end_actual),data(i_ISP+1).time_surf(start_next:end_next)],...
% %                              gap.x_surf(1:N_surf_interp),...
% %                              gap.time_surf(1:N_surf_interp));
% %     
% %     gap.y_surf(:) = spline( [time_surf(start_actual:end_actual),data(i_ISP+1).time_surf(start_next:end_next)],...
% %                              gap.y_surf(1:N_surf_interp),...
% %                              gap.time_surf(1:N_surf_interp));
% %     
% %     gap.z_surf(:) = spline( [time_surf(start_actual:end_actual),data(i_ISP+1).time_surf(start_next:end_next)],...
% %                              gap.z_surf(1:N_surf_interp),...
% %                              gap.time_surf(1:N_surf_interp));
% 
%     %>> Cartesian to geodetic
%     lla = ecef2lla([gap(i_gap).x_sat.',gap(i_gap).y_sat.',gap(i_gap).z_sat.'],flat_coeff_cst,semi_major_axis_cst);
%     gap(i_gap).lat_sat = lla(:,1).';
%     gap(i_gap).lon_sat = lla(:,2).';
%     gap(i_gap).alt_sat = lla(:,3).';
% 
% 
% %>> Orbit to surface
%     gap(i_gap).win_delay_surf = mean([win_delay_surf(end_surf(i_ISP-1)),win_delay_sar_ku(init_burst(i_ISP))])*ones(1,gap(i_gap).N_total_surf_loc);
%     gap(i_gap).lat_surf = gap(i_gap).lat_sat;
%     gap(i_gap).lon_surf = gap(i_gap).lon_sat;
%     gap(i_gap).alt_surf = gap(i_gap).alt_sat - gap(i_gap).win_delay_surf * c_cst/2;
% 
%     %>> Geodetic to cartesian
%     p = lla2ecef([gap(i_gap).lat_surf.',gap(i_gap).lon_surf.',gap(i_gap).alt_surf.'],flat_coeff_cst,semi_major_axis_cst);
%     gap(i_gap).x_surf = p(:,1).';
%     gap(i_gap).y_surf = p(:,2).';
%     gap(i_gap).z_surf = p(:,3).';
% %     gap(i_gap).Process_ID2 = zeros(1,length(gap(i_gap).lon_surf))+min([Process_ID2,data(i_ISP+1).Process_ID2]);
%     
%     
% % %>> Surface to orbit
% %     gap(i_gap).x_sat(:) = spline( [time_sar_ku(start_actual:end_actual),data(i_ISP+1).time_surf(start_next:end_next)],...
% %                              gap(i_gap).x_sat(1:N_surf_interp),...
% %                              gap(i_gap).time_surf(1:N_surf_interp));
% %     
% %     gap(i_gap).y_sat(:) = spline( [time_surf(start_actual:end_actual),data(i_ISP+1).time_surf(start_next:end_next)],...
% %                              gap(i_gap).y_sat(1:N_surf_interp),...
% %                              gap(i_gap).time_surf(1:N_surf_interp));
% %     
% %     gap(i_gap).z_sat(:) = spline( [time_surf(start_actual:end_actual),data(i_ISP+1).time_surf(start_next:end_next)],...
% %                              gap(i_gap).z_sat(1:N_surf_interp),...
% %                              gap(i_gap).time_surf(1:N_surf_interp));
% %     
% %     
% % %>> Cartesian to geodetic
% %     lla = ecef2lla([gap(i_gap).x_sat,gap(i_gap).y_sat,gap(i_gap).z_sat],flat_coeff_cst,semi_major_axis_cst);
% %     gap(i_gap).lat_sat(:) = lla(:,1);
% %     gap(i_gap).lon_sat(:) = lla(:,2);
% %     gap(i_gap).alt_sat(:) = lla(:,3);
% %     
% %     
% % %>> Window delay
% %     gap(i_gap).win_delay_surf = (gap(i_gap).alt_sat - gap(i_gap).alt_surf) * 2 / c_cst;
%     
%     
% %% Setting the first surface location for next ISP
% 
%     init.lat_surf           = gap(i_gap).lat_surf(gap(i_gap).N_total_surf_loc);
%     init.lon_surf           = gap(i_gap).lon_surf(gap(i_gap).N_total_surf_loc);
%     init.alt_surf           = gap(i_gap).alt_surf(gap(i_gap).N_total_surf_loc);
%     init.x_surf             = gap(i_gap).x_surf(gap(i_gap).N_total_surf_loc);
%     init.y_surf             = gap(i_gap).y_surf(gap(i_gap).N_total_surf_loc);
%     init.z_surf             = gap(i_gap).z_surf(gap(i_gap).N_total_surf_loc);
%     init.time_surf          = gap(i_gap).time_surf(gap(i_gap).N_total_surf_loc);
%     init.lat_sat            = gap(i_gap).lat_sat(gap(i_gap).N_total_surf_loc);
%     init.lon_sat            = gap(i_gap).lon_sat(gap(i_gap).N_total_surf_loc);
%     init.alt_sat            = gap(i_gap).alt_sat(gap(i_gap).N_total_surf_loc);
%     init.x_vel_sat          = gap(i_gap).x_vel_sat(gap(i_gap).N_total_surf_loc);
%     init.y_vel_sat          = gap(i_gap).y_vel_sat(gap(i_gap).N_total_surf_loc);
%     init.z_vel_sat          = gap(i_gap).z_vel_sat(gap(i_gap).N_total_surf_loc);
%     init.pitch_surf         = gap(i_gap).pitch_surf(gap(i_gap).N_total_surf_loc);
%     init.roll_surf          = gap(i_gap).roll_surf(gap(i_gap).N_total_surf_loc);
%     init.yaw_surf           = gap(i_gap).yaw_surf(gap(i_gap).N_total_surf_loc);
%     init.win_delay_surf     = gap(i_gap).win_delay_surf(gap(i_gap).N_total_surf_loc);
%     init.alt_rate_sat       = gap(i_gap).alt_rate(gap(i_gap).N_total_surf_loc);
%     % Find from which burst to start computing the surface locations
%     diff_time = time_sar_ku(init_burst(i_ISP):init_burst(i_ISP)+N_bursts_cycle_chd) - init.time_surf;
%     [~,i_curr] = min(abs(diff_time));
%     init.i_curr = i_curr - 1; %we subtract '1' because before we start using it, we add '1' to the value (this is how the loop works best)    
%     
%     %% INSERT bursts info BETWEEN end_burst(i_ISP-1) init_burst(i_ISP) 
%     %% and surf info after end_surf(i_ISP-1)
%     
%     N_total_surf_loc = [N_total_surf_loc gap(i_gap).N_total_surf_loc];
%     
%     
%     
%     
%     
% end
% 
% 
% 
% 
