% compute_height_rate
% It computes an approximation of the height rate (in m/s)

function [height_rate,a] = compute_height_rate(N_total_bursts_sar_ku_isp,...
                                           x_vel_sat_geoloc, y_vel_sat_geoloc, z_vel_sat_geoloc, ...
                                           x_sat_geoloc, y_sat_geoloc, z_sat_geoloc, ...
                                           x_surf_geoloc, y_surf_geoloc, z_surf_geoloc)

global pi_cst
a = zeros(1,N_total_bursts_sar_ku_isp);
height_rate = zeros(1,N_total_bursts_sar_ku_isp);

for i_burst = 1:N_total_bursts_sar_ku_isp
    % nadir direction
    n = [x_surf_geoloc(i_burst) - x_sat_geoloc(i_burst), ...
         y_surf_geoloc(i_burst) - y_sat_geoloc(i_burst), ...
         z_surf_geoloc(i_burst) - z_sat_geoloc(i_burst)];
    n_norm = n/norm(n);
    
    % velocity vector
    v = [x_vel_sat_geoloc(i_burst),y_vel_sat_geoloc(i_burst),z_vel_sat_geoloc(i_burst)];
    v_norm = v/norm(v);

    % Vector perpendicular a 'n' i 'v': 'w' (vector associat al pla A)
    w = cross(n,v);
    w_norm = w/norm(w);

    % Vector perpendicular a 'n' i 'w': 'm'
    m = cross(w_norm,n_norm);
    m_norm = m/norm(m);

   % angle between 'v' and 'm' 
   if acos(dot(v_norm,n_norm))<pi_cst/2 
       a(i_burst) = -1.0*acos(v_norm*m_norm'); 
   else
       a(i_burst) = 1.0*acos(v_norm*m_norm'); 
   end
   
    % height rate
    height_rate(i_burst) = norm(v) * sin(a(i_burst));
end
