
h= figure; subplot(1,3,1);
zp_fact_trp=1;N_samples=128;k=surf(abs((fft(squeeze(out_FBR_wv.wfm_cal_gain_corrected(721,:,:))'))).'.^2); colormap(jet);
figlabels('Samples','Pulses','','Closest Burst',16);
set(gca,'XLim',[1 128]);
set(gca,'YLim',[1 64]);

subplot(1,3,2);
k=surf(waveform); colormap(hot);
set(k, 'edgecolor','none');
figlabels('Samples','Beams','','TRP Stack',16);
set(gca,'XLim',[1 size(waveform,2)]);
set(gca,'YLim',[1  size(waveform,1)]);

subplot(1,3,3); plot(out_L1B(record).waveform_group(DB).averaged_power_echo_waveform);
figlabels('Samples','Power','','TRP L1B',16);
set(gca,'XLim',[1 size(out_L1B(record).waveform_group(DB).averaged_power_echo_waveform,2)]);