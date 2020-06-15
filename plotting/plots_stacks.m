%%%plots poster stacks
global N_samples_sar_chd zp_fact_range_cnf T0_chd c_cst
% lla2kml('./results/CS_OFFL_SIR_SIN_FR_20140818T053651_20140818T053937_B001.kml', L1BS.lat_surf, L1BS.lon_surf, L1BS.alt_surf,'.')
i_surf=150;
%% L1B-S
beams=find(L1BS.ProcessID_beam(i_surf,:)==58);

% i_surf=100;
index_inside= 176:179
for i_surf=min(index_inside):max(index_inside)
% figure; subplot(2,4,1); k=surf((squeeze(L1BS.beams_rng_cmpr(i_surf,:,:))));
% figlabels('Samples','Beams','','',16);colormap(colormapATDD);colorbar;freezeColors;
% set(k, 'edgecolor','none');view(-45,45);
% set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
% subplot(2,4,2); k=surf((squeeze(L1BS.beams_rng_cmpr_2(i_surf,:,:))));
% figlabels('Samples','Beams','','',16);colormap(colormapATDD);colorbar;freezeColors;
% set(k, 'edgecolor','none');view(-45,45);
% set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
 h=figure;subplot(1,2,1); imagesc((squeeze(L1BS.beams_rng_cmpr(i_surf,:,:))));
    figlabels('Samples','Beams','',['Power Stack A1 ' num2str(i_surf, '%03d')],16);colormap('jet');colorbar;freezeColors;
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);

for i_beam = 1:L1BS.N_beams_stack(i_surf)
aux(i_beam,1:N_samples_sar_chd*zp_fact_range_cnf)  = squeeze(L1BS.beams_rng_cmpr(i_surf,i_beam,:).*(L1BS.stack_mask(i_surf,i_beam,:)-1))./max(squeeze(L1BS.beams_rng_cmpr(i_surf,i_beam,:)));
aux2(i_beam,1:N_samples_sar_chd*zp_fact_range_cnf)  = squeeze(L1BS.beams_rng_cmpr_filt(i_surf,i_beam,:))./max(squeeze(L1BS.beams_rng_cmpr_filt(i_surf,i_beam,:)));
        


end
    
    h=figure;subplot(1,2,1); imagesc(aux(beams,:));
    figlabels('Samples','Beams','',['Power Stack A1 ' num2str(i_surf, '%03d')],16);colormap('jet');colorbar;freezeColors;
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    
    subplot(1,2,2); imagesc(squeeze(L1BS.beams_rng_cmpr_filt(i_surf,:,:)));
    figlabels('Samples','Beams','',['Power Stack A1 ' num2str(i_surf, '%03d')],16);colormap('jet');colorbar;freezeColors;
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
   

    %% POWER
    h=figure; subplot(2,2,1); imagesc((squeeze(L1BS.beams_rng_cmpr(i_surf,beams,:))));
    figlabels('Samples','Beams','',['Power Stack A1 ' num2str(i_surf, '%03d')],16);colormap('jet');colorbar;freezeColors;
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    subplot(2,1,2); plot(L1B.wfm_cor_i2q2_sar_ku(i_surf,:));figlabels('Samples','','','',16);
    figlabels('Samples','','',['Power L1B waveforms ' num2str(i_surf, '%03d')],16);
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
    subplot(2,2,2); imagesc((squeeze(L1BS.beams_rng_cmpr_2(i_surf,beams,:))));
    figlabels('Samples','Beams','',['Power Stack A2 ' num2str(i_surf, '%03d')],16);colormap(colormapATDD);colorbar;freezeColors;
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    subplot(2,1,2); hold all;plot(L1B.wfm_cor_i2q2_sar_ku_2(i_surf,:)); legend('A1','A2');
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);   
    saveas (h, ['./results/Power_stack_' num2str(i_surf, '%04d') '.jpg']);
    close(h);
    h=figure;subplot(2,2,1); imagesc((squeeze(L1BS.coherence(i_surf,:,:)))); 
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    figlabels('Samples','Beams','',['Coherence Stack ' num2str(i_surf, '%03d')],16);colormap('hot');freezeColors;
    subplot(2,2,3);plot(L1B.coherence(i_surf,:)); 
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
    figlabels('Samples','','',['Coherence L1B ' num2str(i_surf, '%03d')],16);
    subplot(2,2,2); imagesc((squeeze(L1BS.phase_diff(i_surf,:,:))));
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    figlabels('Samples','Beams','',['Phase diff Stack ' num2str(i_surf, '%03d')],16);colormap('hsv');colorbar;    
    subplot(2,2,4);plot(L1B.phase_diff(i_surf,:));
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
    figlabels('Samples','','',['Phase L1B ' num2str(i_surf, '%03d')],16);
	saveas (h, ['./results/Coherence_stack_' num2str(i_surf, '%04d') '.jpg']);
    close(h);
%     subplot(2,2,2); imagesc((squeeze(L1BS.coherence_mask(i_surf,:,:))));
%     set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)])
%     figlabels('Samples','Beams','','',16);colormap(colormapMASK);freezeColors;    
%     
    
    
    
end

for i_surf=136:165%min(index_inside):max(index_inside)
    
    h=figure('Position',[1,1,1920*3,1080]); subplot(2,4,1); imagesc((squeeze(L1BS.beams_rng_cmpr(i_surf,beams,:))));
    figlabels('Samples','Beams','',['Power Stack A1 ' num2str(i_surf, '%03d')],16);colormap('jet');freezeColors;
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    subplot(2,2,3); plot(L1B.wfm_cor_i2q2_sar_ku(i_surf,:));figlabels('Samples','','','',16);
    figlabels('Samples','','',['Power L1B waveforms ' num2str(i_surf, '%03d')],16);
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
    subplot(2,4,2); imagesc((squeeze(L1BS.beams_rng_cmpr_2(i_surf,:,:))));
    figlabels('Samples','Beams','',['Power Stack A2 ' num2str(i_surf, '%03d')],16);colormap('hot');freezeColors;
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    subplot(2,2,3); hold all;plot(L1B.wfm_cor_i2q2_sar_ku_2(i_surf,:)); legend('A1','A2');
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]); 
    subplot(2,4,3); imagesc((squeeze(L1BS.coherence(i_surf,:,:)))); 
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    figlabels('Samples','Beams','',['Coherence Stack ' num2str(i_surf, '%03d')],16);colormap('hot');freezeColors;
    subplot(2,4,7);plot(L1B.coherence(i_surf,:)); 
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
    figlabels('Samples','','',['Coherence L1B ' num2str(i_surf, '%03d')],16);
    subplot(2,4,4); imagesc((squeeze(L1BS.phase_diff(i_surf,:,:))));
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
    figlabels('Samples','Beams','',['Phase diff Stack ' num2str(i_surf, '%03d')],16);colormap('hsv');colorbar;    
    subplot(2,4,8);plot(L1B.phase_diff(i_surf,:));
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
    figlabels('Samples','','',['Phase L1B ' num2str(i_surf, '%03d')],16);
    saveas (h, ['./results/Stack_Pack_' num2str(i_surf, '%04d') '.jpg']);
    close(h);

end







figure; subplot(3,3,1); imagesc((squeeze(L1BS.beams_rng_cmpr(i_surf,:,:))));
figlabels('Samples','Beams','','L1BS.beams rng cmpr',16);colormap('hot');freezeColors;
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
subplot(3,3,2); imagesc((squeeze(L1BS.beams_rng_cmpr_2(i_surf,:,:))));
figlabels('Samples','Beams','','L1BS.beams rng cmpr 2',16);colormap('hot');freezeColors;
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
subplot(3,3,3); imagesc((squeeze(L1BS.coherence_mask(i_surf,:,:))));
figlabels('Samples','Beams','','L1BS.coherence mask',16);colormap(colormapMASK);freezeColors;
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
subplot(3,3,4); imagesc((squeeze(L1BS.beams_rng_cmpr_Coh(i_surf,:,:))));
figlabels('Samples','Beams','','L1BS.beams rng cmpr Coh',16);colormap('hot');freezeColors;
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);
subplot(3,3,5); imagesc((squeeze(L1BS.beams_rng_cmpr_2_Coh(i_surf,:,:))));
figlabels('Samples','Beams','','L1BS.beams rng cmpr 2 Coh',16);colormap('hot');freezeColors;
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)]);

subplot(3,3,7); plot(L1B.wfm_cor_i2q2_sar_ku(i_surf,:));
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
hold all; plot(L1B.wfm_cor_i2q2_sar_ku_Coh(i_surf,:));figlabels('Samples','','','wfm vs wfm Coh',16);
subplot(3,3,8); plot(L1B.wfm_cor_i2q2_sar_ku_2(i_surf,:));
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
hold all; plot(L1B.wfm_cor_i2q2_sar_ku_2_Coh(i_surf,:));figlabels('Samples','','','wfm vs wfm Coh',16);


%% PHASE
figure;subplot(2,1,1); imagesc((squeeze(L1BS.phase_diff(i_surf,:,:))));
figlabels('Samples','Beams','','',16);colormap('jet');colorbar; freezeColors;
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)])
subplot(2,1,2); plot(L1B.phase_diff(i_surf,:)); 
hold all; plot(L1B.phase_diff_2(i_surf,:));
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
figlabels('Samples','','','',16);

%% COHERENCE
figure;subplot(2,2,1); imagesc((squeeze(L1BS.coherence_mask(i_surf,:,:))));colormap(colormapMASK);freezeColors;
subplot(2,2,2); imagesc((squeeze(L1BS.coherence(i_surf,:,:))));
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 L1BS.N_beams_stack(i_surf)])
figlabels('Samples','Beams','','',16);(colormap('hot'));colorbar;freezeColors;

subplot(2,3,4);plot(L1B.coherence(i_surf,:));
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
figlabels('Samples','','','L1B.coherence',16);
subplot(2,3,5);plot(L1B.coherence_2(i_surf,:));
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
figlabels('Samples','','','L1B.coherence_2',16);
subplot(2,3,6);plot(L1B.coherence_Coh(i_surf,:));
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);
figlabels('Samples','','','L1B.coherence_Coh',16);
% subplot(2,2,1);colormap('hot');




figure;subplot(2,1,2);plot(sum(squeeze(L1BS.phase_diff(i_surf,:,:).*L1BS.coherence_mask(i_surf,:,:))));colormap('hot')
subplot(2,1,2);plot(sum(squeeze(L1BS.phase_diff(i_surf,:,:).*L1BS.coherence_mask(i_surf,:,:))));colormap('hot')



