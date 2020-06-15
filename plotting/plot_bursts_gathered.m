h=figure;

for i_burst = L1BS_buffer(i_surf_stacked).burst_index(1)+40:L1BS_buffer(i_surf_stacked).burst_index(end)-40

	mesh(abs(fft(((squeeze(L1A_buffer(i_burst).beams_focused_shifted(:,:)))).', N_samples ).').^2); hold all;
    colormap(jet); colorbar;
    beam_selected=(abs(fft(((squeeze(L1A_buffer(i_burst).beams_focused_shifted(L1BS_buffer(i_surf_stacked).beam_index(find(L1BS_buffer(i_surf_stacked).burst_index==i_burst)),:)))).', N_samples).').^2);
    plot3(1:128,zeros(1,128)+ L1BS_buffer(i_surf_stacked).beam_index(find(L1BS_buffer(i_surf_stacked).burst_index==i_burst)),reshape(beam_selected,[size(beam_selected) 1]),'or');
    set(gca,'XLim',[1 128],'FontSize',12);
    set(gca,'YLim',[1 64],'FontSize',12);
    figlabels('Samples','Beams','Power',['Burst ' num2str(i_burst)],12)
    view(60,50)
    hold off;
    saveas (h,['./results/plots/Burst_' num2str(i_burst, '%03d') '.png']);

    
    
    
end
set(0,'defaultFigureVisible','on');
zp_fact_trp=2;
h=figure;
subplot(2,2,1);
% k=surf(1:1/zp_fact_trp:N_samples+1-1/zp_fact_trp,1:L1BS.N_beams_stack,abs(fftshift(ifft(L1BS.beams_surf(:,:).',N_samples*zp_fact_trp),1)).'.^2);
k=surf(1:1/zp_fact_trp:N_samples+1-1/zp_fact_trp,1:L1BS.N_beams_stack,abs(fftshift(fft(L1BS.beams_surf(:,:).',N_samples*zp_fact_trp),1)).'.^2);

% k=surf(1:1/zp_fact_trp:N_samples+1-1/zp_fact_trp,1:L1BS.N_beams_stack,abs(fft(L1BS.beams_surf(:,:).',N_samples*zp_fact_trp)).'.^2);
set(k, 'edgecolor','none');view(10,50);
set(gca,'XLim',[1 N_samples],'FontSize',12);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',12);
figlabels('Samples','Beams','Power',['Stack #' num2str(L1BS.surf_counter) ' before alignment'],12)

subplot(2,2,3);
plot(L1BS.shift,1:L1BS.N_beams_stack);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',12);
set(gca,'XLim',[-N_samples N_samples],'FontSize',12);
figlabels('Shift [samples]','Beams','','Range corrections',12)
subplot(2,2,2);
% k=surf(1:1/zp_fact_trp:N_samples+1-1/zp_fact_trp,1:L1BS.N_beams_stack,abs(fftshift(ifft(L1BS.beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2);
k=surf(1:1/zp_fact_trp:N_samples+1-1/zp_fact_trp,1:L1BS.N_beams_stack,abs(fftshift(fft(L1BS.beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).'.^2);
set(k, 'edgecolor','none');
view(10,50);
set(gca,'XLim',[1 N_samples],'FontSize',12);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',12);
figlabels('Samples','Beams','Power',['Stack #' num2str(L1BS.surf_counter) ' aligned'],12)

subplot(2,2,4);
zp_fact_trp=256;
% [max_val,max_pos]=max(abs(fftshift(ifft(L1BS.beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).^2);
[max_val,max_pos]=max(abs(fftshift(fft(L1BS.beam_geo_corr(:,:).',N_samples*zp_fact_trp),1)).^2);
max_pos=max_pos/zp_fact_trp-1/zp_fact_trp+1;
init_beam=42;
end_beam=193;
end_beam=198;
plot(max_pos,1:L1BS.N_beams_stack);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',12);
set(gca,'XLim',[max_pos(floor(L1BS.N_beams_stack/2))-0.05 max_pos(floor(L1BS.N_beams_stack/2))+0.05],'FontSize',12);
set(gca,'XLim',[1 N_samples],'FontSize',12);

[Slope_coef_A1,SS,MUMU] = polyfit(init_beam:end_beam,max_pos(init_beam:end_beam),1);
Slope_A1 = polyval(Slope_coef_A1, init_beam:end_beam,SS,MUMU); %Slope_coef_A1(1) related with datation error
Stack_noise_A1=std(max_pos(init_beam:end_beam)-Slope_A1); % units [mm]
Stack_alignment_A1=(Slope_A1(end)-Slope_A1(1)); % units [mm]
hold all;
plot(Slope_A1,(init_beam:end_beam))
figlabels('Samples','Beams','',['Stack #' num2str(L1BS.surf_counter) ' alignment: ' num2str(Stack_alignment_A1) ' [m] and noise: ' num2str(Stack_noise_A1) ' [m]; zero padding: ' num2str(zp_fact_trp)],12)
hold off