h=figure;
N_samples=256;
L1BS_fid = netcdf.open('C:\Users\Joe\Documents\Sentinel-3 Transponder\Mark checking\L1BS_measurement_l1a_trp.nc');
%L1BS_fid = netcdf.open('I:\S3_TRP_20170227_input_cnf\inputs\measurement_l1a.nc');

lat_trp = 35337930.2808;
lon_trp = 23779518.2869;

lat_lon_diff = zeros(741,3);

for i_surf = 411:411

%[L1BS] = readanyNETCDF_record_L1BS(L1BS_fid,42308-i_surf+100);
[L1BS] = readanyNETCDF_record_L1BS(L1BS_fid,i_surf);
%[L1A] = 
%{
lat_lon_diff(i_surf,1) = lat_trp - L1BS.data.lat_l1bs_echo_sar_ku;
lat_lon_diff(i_surf,2) = lon_trp - L1BS.data.lon_l1bs_echo_sar_ku;
lat_lon_diff(i_surf,3) = mean(lat_lon_diff(i_surf,1:2));
%}

subplot(2,1,1)
k=surf(0:N_samples-1,1:240,abs((((double(L1BS.data.i_echoes_ku_l1bs_echo_sar_ku)+1i*double(L1BS.data.q_echoes_ku_l1bs_echo_sar_ku)))')).^2);
set(k, 'edgecolor','none');view(10,50);

subplot(2,1,2)
plot(double(L1BS.data.i2q2_meas_ku_l1bs_echo_sar_ku));
saveas (h,['Stacks_' num2str(i_surf) '.png']);

end

beam_geo_corr=ifft(fftshift((double(L1BS.data.i_echoes_ku_l1bs_echo_sar_ku)+1i*double(L1BS.data.q_echoes_ku_l1bs_echo_sar_ku)),1)).';

k=imagesc(abs(fftshift(fft(beam_geo_corr.',N_samples*128),1)).'.^2);
set(k, 'edgecolor','none');view(10,50);