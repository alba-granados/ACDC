% figure; plot(L1BS_1.lat_surf,L1BS_1.alt_surf,'-.');
% hold all; plot(L1BS_2.lat_surf,L1BS_2.alt_surf,'o');
% 
% 
% figure; plot(L1BS_1.lat_surf,L1BS_1.alt_surf,'-.');
% hold all; plot(L1BS_2.lat_surf,L1BS_2.alt_surf,'o');
% 
% figure; plot(L1BS_1.shift(180,:)-L1BS_2.shift(359,:));
% 
% 
% figure; plot(L1BS_1.beam_ang_surf(180,:),'o');
% hold all; plot(L1BS_2.beam_ang_surf(359,:),'.');
% 
% 
% figure; plot(L1BS_1.beam_ang_surf(180,:)-L1BS_2.beam_ang_surf(359,:));
% 
% 



% L1BS_1.burst_index(180,:)-L1BS_2.burst_index(359,:)

figure;
az_zp=1;
rg_zp=1;
z_lim=0;
for i_burst=L1BS_1.burst_index(180,1):max(L1BS_1.burst_index(180,:))
    range_fft=abs(fftshift(fft(squeeze(L1A_1.beams_focused_shifted(i_burst,:,:)).',128*rg_zp),1));
    z_lim = max(range_fft,z_lim);
end
for i_burst=L1BS_1.burst_index(180,1):max(L1BS_1.burst_index(180,:))
    range_fft=abs(fftshift(fft(squeeze(L1A_1.beams_focused_shifted(i_burst,:,:)).',128*rg_zp),1));
%     azimuth_fft=abs(fftshift(fft(range_fft.',64*1),1)); 
    mesh(1/rg_zp:1/rg_zp:128,1/az_zp:1/az_zp:64,range_fft.'); view(60,75);
    set(gca,'XLim',[1 128],'Fontsize', 20);
    set(gca,'YLim',[1 64], 'Fontsize', 20);
    set(gca,'ZLim',[0 max(max(z_lim))], 'Fontsize', 20);
    figlabels('Samples','Beams','Power',['Burst ' num2str(i_burst,'%04d')] ,20);colorbar;
    caxis([0 max(max(z_lim))]);
    figName = ['AZ_ZP_1_Scaled_PT_Burst_',num2str(i_burst,'%04d')];
    saveas (gcf,[figName,'.png'])
end

z_lim=0;

for i_burst=L1BS_2.burst_index(359,1):max(L1BS_2.burst_index(359,:))
    range_fft=abs(fftshift(fft(squeeze(L1A_2.beams_focused_shifted(i_burst,:,:)).',128*rg_zp),1));
    z_lim = max(range_fft,z_lim);
end

az_zp=2;
rg_zp=1;
for i_burst=L1BS_2.burst_index(359,1):max(L1BS_2.burst_index(359,:))
    range_fft=abs(fftshift(fft(squeeze(L1A_2.beams_focused_shifted(i_burst,:,:)).',128*rg_zp),1));
%     azimuth_fft=abs(fftshift(fft(range_fft.',64*1),1)); 
    mesh(1/rg_zp:1/rg_zp:128,1/az_zp:1/az_zp:64,range_fft.'); view(60,75);
    set(gca,'XLim',[1 128],'Fontsize', 20);
    set(gca,'YLim',[1 64], 'Fontsize', 20);    
    set(gca,'ZLim',[0 max(max(z_lim))], 'Fontsize', 20);
    figlabels('Samples','Beams','Power',['Burst ' num2str(i_burst,'%04d')] ,20);colorbar;
    caxis([0 max(max(z_lim))]);
    figName = ['AZ_ZP_2_Scaled_PT_Burst_',num2str(i_burst,'%04d')];
    saveas (gcf,[figName,'.png'])
end