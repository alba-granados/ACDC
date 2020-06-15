[L1A]=readanyNETCDF('C:\Users\albert.ISARDSAT\Documents\WORK\S3 MPC\DATA\TRP\20160409\S3A_SR_1_SRA_A__20160409T200004_20160409T200027_20160525T065201_0022_046_370______LN3_O_ST_001.SEN3\measurement_l1a.nc');
[L1BS]=readanyNETCDF('C:\Users\albert.ISARDSAT\Documents\WORK\S3 MPC\DATA\TRP\20160409\S3A_SR_1_SRA_BS_20160409T200004_20160409T200027_20160525T065201_0022_046_370______LN3_O_ST_001.SEN3\measurement_l1bs.nc');
lat_a=double(L1A.data.lat_l1a_echo_sar_ku).*L1A.attributes.lat_l1a_echo_sar_ku.scale_factor;
lon_a=double(L1A.data.lon_l1a_echo_sar_ku).*L1A.attributes.lon_l1a_echo_sar_ku.scale_factor;
alt_a=double(L1A.data.alt_l1a_echo_sar_ku).*L1A.attributes.alt_l1a_echo_sar_ku.scale_factor+L1A.attributes.alt_l1a_echo_sar_ku.add_offset;
range_a=double(L1A.data.range_ku_l1a_echo_sar_ku).*L1A.attributes.range_ku_l1a_echo_sar_ku.scale_factor+L1A.attributes.range_ku_l1a_echo_sar_ku.add_offset;
lat_bs=double(L1BS.data.lat_l1bs_echo_sar_ku).*L1BS.attributes.lat_l1bs_echo_sar_ku.scale_factor;
lon_bs=double(L1BS.data.lon_l1bs_echo_sar_ku).*L1BS.attributes.lon_l1bs_echo_sar_ku.scale_factor;
alt_bs=double(L1BS.data.alt_l1bs_echo_sar_ku).*L1BS.attributes.alt_l1bs_echo_sar_ku.scale_factor+L1BS.attributes.alt_l1bs_echo_sar_ku.add_offset;
range_bs=double(L1BS.data.range_ku_l1bs_echo_sar_ku).*L1BS.attributes.range_ku_l1bs_echo_sar_ku.scale_factor+L1BS.attributes.range_ku_l1bs_echo_sar_ku.add_offset;

[~,pos]=(max(L1A.data.i2q2_meas_ku_l1a_echo_plrm));
semi_major_axis_cst = 6378137;
flat_coeff_cst = 0.00335281066;
p = lla2ecef([lat_a,lon_a,alt_a-range_a],flat_coeff_cst,semi_major_axis_cst);

lat_trp = 35.3379302808;
lon_trp = 23.7795182869;
alt_trp = 1048.8184;

p_trp = lla2ecef([lat_trp,lon_trp,alt_trp],flat_coeff_cst,semi_major_axis_cst);

for i_surf=1:length(lat_a)
    Dist_a(i_surf) = norm([p(i_surf,1), p(i_surf,2), p(i_surf,3)] - [p_trp(1),p_trp(2),p_trp(3)]);
end

for i_surf=1:length(lat_a)
    Dist_a(i_surf) = norm([lat_a(i_surf), lon_a(i_surf), alt_a(i_surf)-range_a(i_surf)] - [lat_trp,lon_trp,alt_trp]);
end
for i_surf=1:length(lat_a)
    Dist_a(i_surf) = norm([lat_a(i_surf), lon_a(i_surf)] - [lat_trp,lon_trp]);
end
for i_surf=1:length(lat_bs)
    Dist_bs(i_surf) = norm([lat_bs(i_surf), lon_bs(i_surf)] - [lat_trp,lon_trp]);
end

figure; plot3(lon_bs,lat_bs,1:length(lat_bs),'o');hold all; plot3(lon_a,lat_a,1:length(lat_a),'.');plot3(lon_trp,lat_trp,1,'.');

i_surf_stacked=197;
figure;
for i_burst_trp=L1BS.data.burst_start_ind_l1bs_echo_sar_ku(i_surf_stacked):L1BS.data.burst_stop_ind_l1bs_echo_sar_ku(i_surf_stacked)
    wfm_cal_gain_corrected=(squeeze(double(L1A.data.i_meas_ku_l1a_echo_sar_ku(:,:,i_burst_trp))+1i.*double(L1A.data.q_meas_ku_l1a_echo_sar_ku(:,:,i_burst_trp)))');
    surf(abs((fft(((wfm_cal_gain_corrected))'))')); view(60,50);
    figlabels('Samples','Beams','FFT p.u.',['Burst not steered #' num2str(i_burst_trp)],20);
    set(gca,'XLim',[1 128]);set(gca,'YLim',[1 64]);

%     saveas (gcf,['Burst #' num2str(i_burst_trp) '.png']);
%     set(gca,'ZLim',[0 5e5]);
    saveas (gcf,['Unfocused_Burst #' num2str(i_burst_trp) '.png']);
end

% [a,b]=min(abs(lat'-lat_trp));
 [a,b]=min(abs(lon'-lon_trp))
b=197;
stacks=double(L1BS.data.i_echoes_ku_l1bs_echo_sar_ku(:,:,:)) + 1i.*double(L1BS.data.q_echoes_ku_l1bs_echo_sar_ku(:,:,:));
figure; mesh(abs(squeeze(stacks(:,1:212,198)))');