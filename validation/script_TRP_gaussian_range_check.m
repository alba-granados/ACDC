
figure; 
zp_fact_trp=512;subplot(2,2,1); plot(0:1/zp_fact_trp:N_samples-1/zp_fact_trp,sum(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2),'-')
hold all
zp_fact_trp=16;plot(0:1/zp_fact_trp:N_samples-1/zp_fact_trp,sum(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2),'-');
zp_fact_trp=8;plot(0:1/zp_fact_trp:N_samples-1/zp_fact_trp,sum(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2),'-');
zp_fact_trp=4;plot(0:1/zp_fact_trp:N_samples-1/zp_fact_trp,sum(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2),'-');
zp_fact_trp=2;plot(0:1/zp_fact_trp:N_samples-1/zp_fact_trp,sum(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2),'-');
zp_fact_trp=1;plot(0:1/zp_fact_trp:N_samples-1/zp_fact_trp,sum(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2),'-');
legend('ZP = 512' , 'ZP = 16', 'ZP = 8', 'ZP = 4', 'ZP = 2', 'ZP = 1');
set(gca,'XLim',[0 128]);
zp_fact_trp=512;
set(gca,'YLim',[0 max(sum(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2))]);
figlabels('Range bin','FFT power units','','TRP L1B waveform 2016/08/22',18);

L1b=sum(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(valid_beams,:).',N_samples*zp_fact_trp),1)).'.^2);
figure; plot(0:1/zp_fact_trp:N_samples-1/zp_fact_trp,L1b);hold all

[max_val_L1b,max_pos_L1b]=max(L1b);
range_axis = 1:1/zp_fact_trp:N_samples+1-1/zp_fact_trp;
range_bins = max_pos_L1b-150:max_pos_L1b+150;
plot(range_axis(range_bins),L1b(range_bins),'o');
cfun_bins = fit ((range_bins)', L1b(range_bins)', 'gauss1');

 a_gauss                 = cfun_bins.a1;
        stack_bin_centre       = cfun_bins.b1;
        stack_std               = cfun_bins.c1/2;
        stack_width             = 2*sqrt(2*log(2))*cfun_bins.c1;
        % Compute the characterization parameters:
 power_fitted = a_gauss * exp (-(range_bins - stack_bin_centre).^2 /(2*stack_std.^2));  
 hold all
    plot(range_axis(range_bins),power_fitted);    
        stack_skewness= skewness(power_fitted);
        stack_kurtosis= kurtosis(power_fitted)-3;