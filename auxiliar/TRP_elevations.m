% Compute elevations from L1BS_buffer(i_surf) Dedop 
global c_cst bw_ku_chd N_samples_sar_chd zp_fact_range_cnf N_bursts_cycle_chd
global alt_trp lon_trp lat_trp flat_coeff_cst semi_major_axis_cst
zp_fact_range_cnf=128;
i_surf          = i_surf_stacked;
%TRP_int_delay   = 6.530;
TRP_int_delay   = 9.88;
% TRP_int_delay   = 0;

TRP_record = floor(i_surf/N_bursts_cycle_chd);

p = lla2ecef([lat_trp,lon_trp,alt_trp],flat_coeff_cst,semi_major_axis_cst);
    x_TRP = p(:,1).';
    y_TRP = p(:,2).';
    z_TRP = p(:,3).';
TRP_coord = [x_TRP,y_TRP,z_TRP];
TRP_coord = [L1BS_buffer(i_surf).x_surf,L1BS_buffer(i_surf).y_surf,L1BS_buffer(i_surf).z_surf];


[~,beam_index_nadir]=min(abs(-pi/2+L1BS_buffer(i_surf).beam_ang_surf'));

% burst_index_nadir = zeros(1,L1BS_buffer(i_surf).N_total_surf_loc);
% elevation = zeros(1,L1BS_buffer(i_surf).N_total_surf_loc);

    burst_index_nadir= L1BS_buffer(i_surf).burst_index(beam_index_nadir);


% figure; mesh((((abs(fftshift(fft(squeeze(L1BS_buffer(i_surf).beam_geo_corr(:,:)).',N_samples_sar_chd*zp_fact_range_cnf),1)).').^2))); colormap(jet);
% set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf],'Fontsize', 20);
%     set(gca,'YLim',[1 64*4], 'Fontsize', 20);
% %    set(gca,'ZLim',[1 14e8], 'Fontsize', 20);
%     figlabels('Samples','Beams','Power','' ,20);colorbar;
% % figure; mesh(((((squeeze(L1BS_buffer(i_surf).beams_rng_cmpr(:,:)))))))
% figure;  plot(sum((((abs(fftshift(fft(squeeze(L1BS_buffer(i_surf).beam_geo_corr(:,:)).',N_samples_sar_chd*zp_fact_range_cnf),1)).').^2))));
% figlabels('Samples','Power','','' ,20);




geophysical_corrections_sea_ice  = L1A_buffer(burst_index_nadir).dry_tropo_correction_bursts+L1A_buffer(burst_index_nadir).wet_tropo_correction_bursts+ L1A_buffer(burst_index_nadir).inverse_baro_correction_bursts+           L1A_buffer(burst_index_nadir).GIM_iono_correction_bursts+ L1A_buffer(burst_index_nadir).ocean_equilibrium_tide_bursts+  L1A_buffer(burst_index_nadir).long_period_tide_height_bursts+   L1A_buffer(burst_index_nadir).ocean_loading_tide_bursts+ L1A_buffer(burst_index_nadir).solid_earth_tide_bursts+ L1A_buffer(burst_index_nadir).geocentric_polar_tide_bursts;
geophysical_corrections_land_ice = L1A_buffer(burst_index_nadir).dry_tropo_correction_bursts+L1A_buffer(burst_index_nadir).wet_tropo_correction_bursts+                                                                         L1A_buffer(burst_index_nadir).GIM_iono_correction_bursts+                                                                                                                               L1A_buffer(burst_index_nadir).ocean_loading_tide_bursts+ L1A_buffer(burst_index_nadir).solid_earth_tide_bursts+ L1A_buffer(burst_index_nadir).geocentric_polar_tide_bursts;
geophysical_corrections_ocean    = L1A_buffer(burst_index_nadir).dry_tropo_correction_bursts+L1A_buffer(burst_index_nadir).wet_tropo_correction_bursts+ L1A_buffer(burst_index_nadir).Dynamic_atmospheric_correction_bursts+    L1A_buffer(burst_index_nadir).GIM_iono_correction_bursts+ L1A_buffer(burst_index_nadir).ocean_equilibrium_tide_bursts+  L1A_buffer(burst_index_nadir).long_period_tide_height_bursts+   L1A_buffer(burst_index_nadir).ocean_loading_tide_bursts+ L1A_buffer(burst_index_nadir).solid_earth_tide_bursts+ L1A_buffer(burst_index_nadir).geocentric_polar_tide_bursts;

% figure; plot(geophysical_corrections_sea_ice); hold all; plot(geophysical_corrections_land_ice); plot(geophysical_corrections_ocean);
[max_val,max_pos_s]= max(((abs(fftshift(fft(squeeze(L1BS_buffer(i_surf).beam_geo_corr(:,:)).',N_samples_sar_chd*zp_fact_range_cnf),1))).^2));

% figure; plot(L1BS_buffer(i_surf).alt_surf-(((max_pos'./zp_fact_range_cnf-N_samples_sar_chd/2+1)/bw_ku_chd)*c_cst/2) + TRP_int_delay/2-geophysical_corrections_land_ice(floor(burst_index_nadir/20)+1)); % The alt_surf variable is the difference between the satellite altitude and the window delay.
window_corr = L1BS_buffer(i_surf).win_delay_surf-TRP_int_delay/c_cst;

time_stack = window_corr -((((max_pos_s-1)'./zp_fact_range_cnf+1)-(N_samples_sar_chd/2))/bw_ku_chd);

range_stack = time_stack *c_cst/2  + geophysical_corrections_land_ice;
range_without_slant =range_stack-(L1BS_buffer(i_surf).slant_range_corr(:))/bw_ku_chd*c_cst/2;

for i_beam=1:length(range_without_slant)
    
    SAT_coord = [L1BS_buffer(i_surf).x_sar_sat_beam(i_beam),L1BS_buffer(i_surf).y_sar_sat_beam(i_beam),L1BS_buffer(i_surf).z_sar_sat_beam(i_beam)];
    dif       = (SAT_coord(1)-TRP_coord(1))^2 + (SAT_coord(2)-TRP_coord(2))^2 + (SAT_coord(3)-TRP_coord(3))^2;
    theo_range(i_beam) =  sqrt(dif);
end
figure; subplot(3,2,1);
imagesc((((abs(fftshift(fft(squeeze(L1BS_buffer(i_surf).beam_geo_corr(:,:)).',N_samples_sar_chd*zp_fact_range_cnf),1)).').^2))); colormap(jet);
figlabels('Samples','Beams','Power','TRP L1BS' ,14);
subplot(3,2,3);
plot(sum((((abs(fftshift(fft(squeeze(L1BS_buffer(i_surf).beam_geo_corr(:,:)).',N_samples_sar_chd*zp_fact_range_cnf),1)).').^2))));
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
figlabels('Samples','Power','','TRP L1B' ,14);
subplot(3,2,5);
plot(max_pos_s,1:length(range_without_slant))
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
figlabels('Beams','Sample index','','TRP peak position' ,14);
subplot(2,2,2);
plot(range_without_slant'); hold all; plot(theo_range');
figlabels('Beams','Range [m]','','Range' ,14);
legend('Range meas', 'Range theo');
subplot(2,2,4);plot(range_without_slant'-theo_range);
set(gca,'YLim',[-0.4 0.4]);
figlabels('Beams','Range [m]','','Error = Range meas - Range theo' ,14);



[max_val,max_pos]= max(sum((((abs(fftshift(fft(squeeze(L1BS_buffer(i_surf).beam_geo_corr(:,:)).',N_samples_sar_chd*zp_fact_range_cnf),1)).').^2))));
subplot(3,2,1);
set(gca,'XLim',[max_pos-3*zp_fact_range_cnf max_pos+3*zp_fact_range_cnf]);
subplot(3,2,3);
set(gca,'XLim',[max_pos-3*zp_fact_range_cnf max_pos+3*zp_fact_range_cnf]);
subplot(3,2,5);
set(gca,'XLim',[max_pos-3*zp_fact_range_cnf max_pos+3*zp_fact_range_cnf]);

uncorr_elevation = L1BS_buffer(i_surf).alt_surf-(((((max_pos-1)'./zp_fact_range_cnf+1)-(N_samples_sar_chd/2))/bw_ku_chd)*c_cst/2);
elevation = L1BS_buffer(i_surf).alt_sat-L1BS_buffer(i_surf).win_delay_surf*c_cst/2-(((((max_pos-1)'./zp_fact_range_cnf+1)-(N_samples_sar_chd/2))/bw_ku_chd)*c_cst/2) + TRP_int_delay/2+ geophysical_corrections_land_ice;
   
disp(elevation-alt_trp);