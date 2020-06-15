% Compute elevations from L1B Dedop 
  

%plots waveforms L1B
inputFile='./results/data/CR2_SR_1_SRA____20150102T193036_20150102T193250_20160411_isd.nc';
% inputFile='S6_C1B_SAR_L1B-S_NoAmb.nc';
S6_L1B = readanyNETCDF_V1([inputFile]);
set_default_plot;

bw_ku_chd = 320e6; %2.531645569620253e-09;
c_cst  = 299792458;
time = double(S6_L1B.data.time_l1b_echo_sar_ku);
N_total_surf_loc = length(find(time<max(time)));

zp_fact_range_cnf = 1;
N_samples_sar_chd = length(S6_L1B.data.i2q2_meas_ku_l1b_echo_sar_ku(:,1))./zp_fact_range_cnf;

for i_surf=1:N_total_surf_loc    
    wfm_cor_i2q2_sar_ku_netcdf(i_surf,:) = double(S6_L1B.data.i2q2_meas_ku_l1b_echo_sar_ku(:,i_surf)).*(double(S6_L1B.data.waveform_scale_factor_l1b_echo_sar_ku(i_surf)));
%     wfm_cor_i2q2_sar_ku_netcdf(i_surf,:) = wfm_cor_i2q2_sar_ku_netcdf(i_surf,:)./max(wfm_cor_i2q2_sar_ku_netcdf(i_surf,:));
end
% altimeter_range_calibrated_ku = uint32(((win_delay_surf) * c_cst/2 - 1.3e6) * 1e4);
% com_altitude_ku = uint32((alt_sat - 1.3e6) * 1e4);

win_delay_surf_netcdf = (double(S6_L1B.data.range_ku_l1b_echo_sar_ku.').*double(S6_L1B.attributes.range_ku_l1b_echo_sar_ku.scale_factor)+double(S6_L1B.attributes.range_ku_l1b_echo_sar_ku.add_offset))/c_cst*2;
alt_sat_netcdf = double(S6_L1B.data.alt_l1b_echo_sar_ku.').*double(S6_L1B.attributes.alt_l1b_echo_sar_ku.scale_factor)+double(S6_L1B.attributes.alt_l1b_echo_sar_ku.add_offset);

wd_shift_netcdf = zeros(1,N_total_surf_loc);
wfm_cor_i2q2_sar_ku_wdcorr_netcdf = zeros(N_total_surf_loc,N_samples_sar_chd*4);
dry_tropo_correction    = double (S6_L1B.data.dry_tropo_correction_l1b_echo_sar_ku).*double(S6_L1B.attributes.dry_tropo_correction_l1b_echo_sar_ku.scale_factor);
wet_tropo_correction    = double (S6_L1B.data.wet_tropo_correction_l1b_echo_sar_ku).*double(S6_L1B.attributes.wet_tropo_correction_l1b_echo_sar_ku.scale_factor);
inverse_baro_correction = double (S6_L1B.data.inverse_baro_correction_l1b_echo_sar_ku).*double(S6_L1B.attributes.inverse_baro_correction_l1b_echo_sar_ku.scale_factor);
Dynamic_atmospheric_correction = double (S6_L1B.data.Dynamic_atmospheric_correction_l1b_echo_sar_ku).*double(S6_L1B.attributes.Dynamic_atmospheric_correction_l1b_echo_sar_ku.scale_factor);
GIM_iono_correction = double (S6_L1B.data.GIM_iono_correction_l1b_echo_sar_ku).*double(S6_L1B.attributes.GIM_iono_correction_l1b_echo_sar_ku.scale_factor);
ocean_equilibrium_tide = double (S6_L1B.data.ocean_equilibrium_tide_l1b_echo_sar_ku).*double(S6_L1B.attributes.ocean_equilibrium_tide_l1b_echo_sar_ku.scale_factor);
long_period_tide_height = double (S6_L1B.data.long_period_tide_height_l1b_echo_sar_ku).*double(S6_L1B.attributes.long_period_tide_height_l1b_echo_sar_ku.scale_factor);
ocean_loading_tide = double (S6_L1B.data.ocean_loading_tide_l1b_echo_sar_ku).*double(S6_L1B.attributes.ocean_loading_tide_l1b_echo_sar_ku.scale_factor);
solid_earth_tide = double (S6_L1B.data.solid_earth_tide_l1b_echo_sar_ku).*double(S6_L1B.attributes.solid_earth_tide_l1b_echo_sar_ku.scale_factor);
geocentric_polar_tide = double (S6_L1B.data.geocentric_polar_tide_l1b_echo_sar_ku).*double(S6_L1B.attributes.geocentric_polar_tide_l1b_echo_sar_ku.scale_factor);



burst_index_nadir = zeros(1,N_total_surf_loc);
elevation = zeros(1,N_total_surf_loc);

geophysical_corrections_sea_ice  = dry_tropo_correction+wet_tropo_correction+                                       inverse_baro_correction+    GIM_iono_correction+ ocean_equilibrium_tide+ long_period_tide_height+ ocean_loading_tide+ solid_earth_tide+ geocentric_polar_tide;
geophysical_corrections_land_ice = dry_tropo_correction+wet_tropo_correction+                                                                   GIM_iono_correction+                                                  ocean_loading_tide+ solid_earth_tide+ geocentric_polar_tide;
geophysical_corrections_ocean    = dry_tropo_correction+wet_tropo_correction+ Dynamic_atmospheric_correction+                                   GIM_iono_correction+ ocean_equilibrium_tide+ long_period_tide_height+ ocean_loading_tide+ solid_earth_tide+ geocentric_polar_tide;


for i_surf=1:N_total_surf_loc
    
   [max_val,max_pos(i_surf)]= max(wfm_cor_i2q2_sar_ku_netcdf(i_surf,:));
   uncorr_elevation(i_surf) =alt_sat_netcdf(i_surf)-(((max_pos(i_surf)'/zp_fact_range_cnf- N_samples_sar_chd/2+1)/bw_ku_chd)*c_cst/2);
    elevation(i_surf) = alt_sat_netcdf(i_surf) - win_delay_surf_netcdf(i_surf)*c_cst/2 - (((max_pos(i_surf)'/zp_fact_range_cnf- N_samples_sar_chd/2+1)/bw_ku_chd)*c_cst/2) - geophysical_corrections_land_ice(i_surf);
       
end
 figure;plot(elevation(elevation<1e6))