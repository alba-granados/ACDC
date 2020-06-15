set(0,'defaultFigureVisible','on');
aux = L1A.wfm_cal_gain_corrected;
figure;
subplot(2,2,1); mesh(abs(fftshift(fft(aux.'))));view(90,0);
figlabels('Range Bin','Range Bin','FFT units','L1A burst',12);
subplot(2,2,2); mesh(abs(fftshift(fft(aux.'))));view(0,0);
figlabels('Pulse Index','Range Bin','FFT units','L1A burst',12);
subplot(2,2,3);plot(sum(abs(fftshift(fft(aux.')).')));
figlabels('Range Bin','FFT units','','L1A burst averaged',12);
subplot(2,2,4); plot(sum(abs(fftshift(fft(aux.')))));
figlabels('Pulse Index','FFT units','','L1A burst averaged',12);


aux= L1A_buffer(1).beams_focused_shifted;
figure;
subplot(2,2,1); mesh(abs(fftshift(fft(aux.'))));view(90,0);
figlabels('Range Bin','Range Bin','FFT units','L1A burst',12);
subplot(2,2,2); mesh(abs(fftshift(fft(aux.'))));view(0,0);
figlabels('Pulse Index','Range Bin','FFT units','L1A burst',12);
subplot(2,2,3);plot(sum(abs(fftshift(fft(aux.')).')));
figlabels('Range Bin','FFT units','','L1A burst averaged',12);
subplot(2,2,4); plot(sum(abs(fftshift(fft(aux.')))));
figlabels('Pulse Index','FFT units','','L1A burst averaged',12);

aux= L1BS_buffer(1).beams_surf;
figure;
subplot(2,2,1); mesh(abs(fftshift(fft(aux.'))));view(90,0);
figlabels('Range Bin','Range Bin','FFT units','L1BS stack',12);
subplot(2,2,2); mesh(abs(fftshift(fft(aux.'))));view(0,0);
figlabels('Beam Index','Range Bin','FFT units','L1BS stack',12);
subplot(2,2,3);plot(sum(abs(fftshift(fft(aux.')).')));
figlabels('Range Bin','FFT units','','L1BS stack averaged',12);
subplot(2,2,4); plot(sum(abs(fftshift(fft(aux.')))));
figlabels('Beam Index','FFT units','','L1BS stack averaged',12);


aux= L1BS_buffer(1).beam_geo_corr;
figure;
subplot(2,2,1); mesh(abs(fftshift(fft(aux.'))));view(90,0);
figlabels('Range Bin','FFT units','','L1BS stack',12);
subplot(2,2,2); mesh(abs(fftshift(fft(aux.'))));view(0,0);
figlabels('Beam Index','Range Bin','','L1BS stack',12);
subplot(2,2,3);plot(sum(abs(fftshift(fft(aux.')).')));
figlabels('Range Bin','FFT units','','L1BS stack averaged',12);
subplot(2,2,4); plot(sum(abs(fftshift(fft(aux.')))));
figlabels('Beam Index','FFT units','','L1BS stack averaged',12);

aux= L1BS_buffer(1).beams_rng_cmpr;
figure;
subplot(2,2,1); mesh(abs(((aux.'))));view(90,0);
figlabels('Range Bin','FFT units','','L1BS stack after range compression',12);
subplot(2,2,2); mesh(abs(((aux.'))));view(0,0);
figlabels('Beam Index','Range Bin','','L1BS stack after range compression',12);
subplot(2,2,3);plot(sum(abs(((aux.')).')));
figlabels('Range Bin','FFT units','','L1BS stack  after range compression averaged',12);
subplot(2,2,4); plot(sum(abs(((aux.')))));
figlabels('Beam Index','FFT units','','L1BS stack after range compression averaged',12);
