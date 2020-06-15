% Compute elevations from L1B Dedop 
global c_cst bw_ku_chd N_samples_sar_chd zp_fact_range_cnf 

[~,beam_index_nadir]=min(abs(-pi/2+L1BS.beam_ang_surf'));

burst_index_nadir = zeros(1,L1BS.N_total_surf_loc);
elevation_i = zeros(1,L1BS.N_total_surf_loc);
for i_surf = 1:L1BS.N_total_surf_loc
    burst_index_nadir(i_surf)= L1BS.burst_index(i_surf,beam_index_nadir(i_surf));
end
geophysical_corrections_sea_ice  = L1A.dry_tropo_correction+L1A.wet_tropo_correction+                                       L1A.inverse_baro_correction+    L1A.GIM_iono_correction+ L1A.ocean_equilibrium_tide+ L1A.long_period_tide_height+ L1A.ocean_loading_tide+ L1A.solid_earth_tide+ L1A.geocentric_polar_tide;
geophysical_corrections_land_ice_i = L1A.dry_tropo_correction+L1A.wet_tropo_correction+                                                                       L1A.GIM_iono_correction+                                                          L1A.ocean_loading_tide+ L1A.solid_earth_tide+ L1A.geocentric_polar_tide;
geophysical_corrections_ocean    = L1A.dry_tropo_correction+L1A.wet_tropo_correction+ L1A.Dynamic_atmospheric_correction+                                   L1A.GIM_iono_correction+ L1A.ocean_equilibrium_tide+ L1A.long_period_tide_height+ L1A.ocean_loading_tide+ L1A.solid_earth_tide+ L1A.geocentric_polar_tide;


for i_surf=1:L1BS.N_total_surf_loc-3
    
   [max_val,max_pos_i(i_surf)]= max(L1B.wfm_cor_i2q2_sar_ku(i_surf,:));
   uncorr_elevation(i_surf) = L1BS.alt_surf(i_surf)-(((max_pos_i(i_surf)'/zp_fact_range_cnf-N_samples_sar_chd/2+1)/bw_ku_chd)*c_cst/2);
    elevation_i(i_surf) = L1BS.alt_surf(i_surf)-(((max_pos_i(i_surf)'/zp_fact_range_cnf-N_samples_sar_chd/2+1)/bw_ku_chd)*c_cst/2) - geophysical_corrections_land_ice_i(floor(burst_index_nadir(i_surf)/20)+1);
    
    
end
figure;plot(elevation_i)