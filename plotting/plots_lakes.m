global c_cst N_samples_sar_chd T0_chd zp_fact_range_cnf

Window_size=length(abs(squeeze(L1BS.beam_geo_corr_2(1,:,:)).'));
[lakeTiff, R] = geotiffread(['./inputs/mask/lake2.tif']);
AspectRatio   = R.RasterSize(2)/R.RasterSize(1);
Reduced_Size  = round(R.RasterSize/1.5);
latMargin    = Reduced_Size/2*R.DeltaLat;
lonMargin    = Reduced_Size/2*R.DeltaLon;
% h=figure;subplot(2,3,1); mapshow(lakeTiff, R);hold on;axis image off; axis image on; %set(gca,'XLim',[80.75 81.2]);set(gca,'YLim',[34.4 34.7]);
% [~,first]= min(abs(dataL1_single.lon-L1BS.lon_surf(100)));
% [~,first]= min(abs(dataL1_single.lat-L1BS.lat_surf(1)));

for i_surf = 1:L1BS.N_total_surf_loc
    [~,index(i_surf)] = min(abs(dataL1_single.lat-L1BS.lat_surf(i_surf)));  
end

% for i_surf = 101:L1BS.N_total_surf_loc
%    
%    beams=find(L1BS.ProcessID_beam(i_surf,:)==58);
%    subplot(2,3,2:3);imagesc((abs(squeeze(L1BS.beam_geo_corr(i_surf,beams,:)))));set(gca,'XLim',[1 Window_size]);
%    subplot(2,3,5:6);plot(1:Window_size,sum(abs(squeeze(L1BS.beam_geo_corr(i_surf,beams,:)))));set(gca,'XLim',[1 Window_size]);
%    saveas (h, ['./results/wvf_',num2str(i_surf),'_Lake.jpg']); 
% end
for i_surf = 120:L1BS.N_total_surf_loc 
    h=figure;subplot(2,3,1); mapshow(lakeTiff, R);hold on;axis image off; axis image on;
    subplot(2,3,1);hold on;
    scatter(L1BS.lon_surf(1:i_surf),L1BS.lat_surf(1:i_surf),40,[1 1 0],'fill','LineWidth',1);
    scatter(L1BS.lon_surf(i_surf),L1BS.lat_surf(i_surf),40,[0 1 0],'fill','LineWidth',1);
    set(gca,'YLim',[L1BS.lat_surf(i_surf)+latMargin(1) L1BS.lat_surf(i_surf)-latMargin(1)]);
    set(gca,'XLim',[L1BS.lon_surf(i_surf)-lonMargin(1) L1BS.lon_surf(i_surf)+lonMargin(1)]);
     title('L1B record position');
%     ISF=squeeze(isfinite(L1BS.beam_geo_corr(i_surf,:,:)));
    subplot(2,3,1);
    [win_delay_ref(i_surf),win_delay_ref_index(i_surf)] = min(L1A.win_delay_sar_ku(L1BS.burst_index(i_surf,1:L1BS.N_beams_stack(i_surf))));
        
    first_sample    = L1BS.alt_sat(win_delay_ref_index(i_surf))-(win_delay_ref(i_surf))*c_cst/2+N_samples_sar_chd/2*T0_chd*c_cst/2;
    first_sampleL1  = dataL1_single.alt(index(i_surf))*1e-3-(dataL1_single.winDel(index(i_surf))*1e-12)*c_cst/2+N_samples_sar_chd/2*T0_chd/2*c_cst/2;
    end_sample      = L1BS.alt_sat(win_delay_ref_index(i_surf))-(win_delay_ref(i_surf))*c_cst/2-(Window_size-N_samples_sar_chd/2)*T0_chd*c_cst/2;
    end_sampleL1    = dataL1_single.alt(index(i_surf))*1e-3-(dataL1_single.winDel(index(i_surf))*1e-12)*c_cst/2-(N_samples_sar_chd/2)*T0_chd/2*c_cst/2;
    x_axis = first_sample-T0_chd/zp_fact_range_cnf*c_cst/2:-T0_chd/zp_fact_range_cnf*c_cst/2:end_sample;
    x_axisL1 = first_sampleL1-T0_chd/2*c_cst/2:-T0_chd/2*c_cst/2:end_sampleL1;
    y_axis = L1A.lat_sar_sat(L1BS.burst_index(i_surf,1:L1BS.N_beams_stack(i_surf)));
    subplot(2,3,4);
    plot(x_axisL1,dataL1_single.Wvf(index(i_surf),:)/max(dataL1_single.Wvf(index(i_surf),:))); set(gca,'XLim',[end_sampleL1 first_sampleL1]);set(gca,'XDir','reverse'); 
    figlabels('Elevation [m]','Normalised Power','','L1B Waveform from CR2',12);
     
    beams=find(L1BS.ProcessID_beam(i_surf,:)==58);
    subplot(2,3,2:3);
%     k=surf(y_axis(beams),x_axis,abs(squeeze(L1BS.beam_geo_corr_2(i_surf,beams,:)))');set(k, 'edgecolor','none');view(-90,90);
    imagesc(x_axis,y_axis,(abs(squeeze(L1BS.beam_geo_corr_2(i_surf,beams,:)))));set(gca,'XLim',[end_sample first_sample]);set(gca,'YDir','normal');set(gca,'XDir','reverse');
    figlabels('Elevation [m]','Latitude [degrees]','','Extended Stack data from isardSAT processor',12);
    
    subplot(2,3,5:6);
    plot(x_axis,sum(abs(squeeze(L1BS.beam_geo_corr(i_surf,beams,:))))/max(sum(abs(squeeze(L1BS.beam_geo_corr(i_surf,beams,:))))));set(gca,'XLim',[end_sample first_sample]);set(gca,'XDir','reverse'); 
    figlabels('Elevation [m]','Normalised Power','','Extended L1B Waveform from isardSAT processor',12);
    saveas (h, ['./results/CR2 L2 vs L1 vs L1BS/',dataL1_single.date,'/hamming/',dataL1_single.date,'_',num2str(i_surf, '%03d'),'_L1vsL1mod.jpg']);

%     saveas (h, ['./results/CR2 L2 vs L1 vs L1BS/20140721/',num2str(i_surf, '%03d'),'_StackLake.jpg']);
    
    close (h);
end
    
    