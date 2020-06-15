function plot_stacks_analysis_coherence(L1BS,L1B)

global c_cst T0_chd wv_length_ku
global zp_fact_range_cnf N_samples N_max_beams_stack_chd 
global product_coord


colormapMASK=ones(64,3);
colormapMASK(1,:)=[1 0 0];
colormapMASK(end,:)=[0 1 0];



elev_init=-L1BS.win_delay_sar_ku_beam(L1BS.beam_ref2)*c_cst/2+L1BS.alt_sat_beam(L1BS.beam_ref2);
elev_wd=-L1BS.win_delay_sar_ku_beam(L1BS.beam_ref)*c_cst/2+L1BS.alt_sat_beam(L1BS.beam_ref); 

y_axis = (elev_wd - T0_chd*c_cst/2/zp_fact_range_cnf + N_samples/2*T0_chd*c_cst/2):-T0_chd*c_cst/2/zp_fact_range_cnf:elev_wd-N_samples*T0_chd*c_cst/2+ N_samples/2*T0_chd*c_cst/2;
y_axis_extend = (elev_init - T0_chd*c_cst/2/zp_fact_range_cnf + N_samples/2*T0_chd*c_cst/2):-T0_chd*c_cst/2/zp_fact_range_cnf:elev_init-N_samples*L1BS.N_windows*T0_chd*c_cst/2+ N_samples/2*T0_chd*c_cst/2;

Xdist = L1BS.phase_diff.* wv_length_ku * L1BS.win_delay_surf *c_cst/2./(2*pi*1.1676)/180*pi;
phase_diff_filt=L1BS.phase_diff.*L1BS.coherence_mask;
phase_diff=L1BS.phase_diff;
phase_diff(isnan(L1BS.phase_diff))=0;
phase_diff_filt(isnan(phase_diff_filt))=0;
Xdist_filt=Xdist.*L1BS.coherence_mask;
% Xdist(isnan(Xdist))=0;
Xdist_filt(isnan(Xdist_filt))=0;

% h=figure; 
% colormapSTACKS=colormap(jet(1024*8));
% colormapSTACKS(1,:)=[1 1 1];
% 
% subplot(3,3,1:3); 
%                 scatter(product_coord(1,:),product_coord(4,:),83,'MarkerEdgeColor',[0 0 0]);
%                 hold on; scatter(L1BS.lat_sat_beam,(-L1BS.wd_corr(1:length(L1BS.lat_sat_beam))+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,83,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);%view(10,90); 
%                 hold on; scatter(L1BS.lat_sat_beam(L1BS.beam_ref),(-L1BS.wd_corr(L1BS.beam_ref)+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,80,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','red');
%                 set(gca,'XLim',[min(product_coord(1,product_coord(1,:)~=0)) max(product_coord(1,:))]);
% %                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                 figlabels('Latitude [degrees]','Elevation [m]','','Tracking window position' ,12);
% for i_window=0:L1BS.N_windows
%     hold on; plot(L1BS.lat_sat_beam, zeros(1,length(L1BS.lat_sat_beam))+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
% end
% % coherence_plot=L1BS.coherence_mask;
% % coherence_plot(isnan(coherence_plot))=0;
% % subplot(3,3,[5 8]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),coherence_plot'); colormap(colormapMASK);freezeColors;
% %                         figlabels('Beams','Elevation [m]','','Coherence Mask',12);%view(10,26);
% %                         set(gca,'YDir','normal');
% % %                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)])
% % %                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
% 
% subplot(3,2,[3 5]); k=surf(1:N_max_beams_stack_chd,y_axis_extend,L1BS.beams_rng_cmpr');view(10,26);
%                         set(k, 'edgecolor','none');colormap(colormapSTACKS); freezeColors;
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
%                         figlabels('Beams','Elevation [m]','', ['Stack # ' num2str(L1BS.surf_counter, '%03d') ' before Coherence masking'] ,12);
%                         for i_window=0:L1BS.N_windows
%                             hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%                         end
%         
% subplot(3,2,[4 6]); k=surf(1:N_max_beams_stack_chd,y_axis_extend,L1BS.beams_rng_cmpr_Coh');view(10,26);
%                         set(k, 'edgecolor','none');colormap(colormapSTACKS); freezeColors;
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
%                         figlabels('Beams','Elevation [m]','',['Stack # ' num2str(L1BS.surf_counter, '%03d') ' after Coherence masking'] ,12);
%                         for i_window=0:L1BS.N_windows
%                             hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%                         end
% %                         
% % subplot(2,2,4); plot(y_axis_extend,L1B.wfm_cor_i2q2_sar_ku_Coh,'g');hold all; plot(y_axis_extend,L1B.wfm_cor_i2q2_sar_ku,'b');
% %  figlabels('Elevation [m]','FFT power units','',['L1B waveform # ' num2str(L1BS.surf_counter, '%03d')] ,12);
% %  set(gca,'XLim',[min(y_axis_extend) max(y_axis_extend)]);
% %  legend('Nominal', 'Coherence Filtered')
% saveas (h,['Stack_Extended_Filtered_' num2str(L1BS.surf_counter, '%03d') '.png']);
% close(h);
% 
% %% PHASE 
% h=figure; 
% colormapPHASE = colormap(jet);
% % colormapPHASE(1024/2,:)=[1 1 1];
% % colormapPHASE(1024/2-1,:)=[1 1 1];
% % colormapPHASE(1024/2+1,:)=[1 1 1];
% 
% subplot(3,3,1:3); 
%                 scatter(product_coord(1,:),product_coord(4,:),83,'MarkerEdgeColor',[0 0 0]);
%                 hold on; scatter(L1BS.lat_sat_beam,(-L1BS.wd_corr(1:length(L1BS.lat_sat_beam))+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,83,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);%view(10,90); 
%                 hold on; scatter(L1BS.lat_sat_beam(L1BS.beam_ref),(-L1BS.wd_corr(L1BS.beam_ref)+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,80,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','red');
%                 set(gca,'XLim',[min(product_coord(1,product_coord(1,:)~=0)) max(product_coord(1,:))]);
% %                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                 figlabels('Latitude [degrees]','Elevation [m]','','Tracking window position' ,12);
% for i_window=0:L1BS.N_windows
%     hold on; plot(L1BS.lat_sat_beam, zeros(1,length(L1BS.lat_sat_beam))+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
% end
% 
% % subplot(3,3,[5 8]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),coherence_plot'); colormap(colormapMASK);freezeColors;
% %                         figlabels('Beams','Elevation [m]','','Coherence Mask',12);%view(10,26);
% %                         set(gca,'YDir','normal');
% % %                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)])
% % %                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
% 
% subplot(3,2,[3 5]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),phase_diff');
%                         colormap(colormapPHASE); freezeColors; set(gca,'YDir','normal');
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
%                         figlabels('Beams','Elevation [m]','', ['Phase Difference # ' num2str(L1BS.surf_counter, '%03d') ' before Coherence masking'] ,12);
%                         for i_window=0:L1BS.N_windows
%                             hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%                         end
%                         colorbar;
% 
% subplot(3,2,[4 6]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),phase_diff_filt');
%                         colormap(colormapPHASE); freezeColors; set(gca,'YDir','normal');
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
%                         figlabels('Beams','Elevation [m]','',['Phase Difference # ' num2str(L1BS.surf_counter, '%03d') ' after Coherence masking'] ,12);
%                         for i_window=0:L1BS.N_windows
%                             hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%                         end
%                         colorbar;
% saveas (h,['PhaseDiff_Extended_Filtered_' num2str(L1BS.surf_counter, '%03d') '.png']);
%       
% close(h);
% % 
% % h= figure;  plot(y_axis_extend,L1B.wfm_cor_i2q2_sar_ku_Coh,'g');hold all; plot(y_axis_extend,L1B.wfm_cor_i2q2_sar_ku,'b');
% %  figlabels('Elevation [m]','FFT power units','',['L1B waveform # ' num2str(L1BS.surf_counter, '%03d')] ,12);
% %  set(gca,'XLim',[min(y_axis_extend) max(y_axis_extend)]);
% %  legend('Nominal', 'Coherence Filtered')
% % saveas (h,['L1B_Comparison_' num2str(L1BS.surf_counter, '%03d') '.png']);
% % saveas (h,['L1B_Comparison_' num2str(L1BS.surf_counter, '%03d') '.fig']);
% % 
% % close(h);
% 
% h= figure;
% 
% 
% subplot(3,3,1:3); 
%                 scatter(product_coord(1,:),product_coord(4,:),83,'MarkerEdgeColor',[0 0 0]);
%                 hold on; scatter(L1BS.lat_sat_beam,(-L1BS.wd_corr(1:length(L1BS.lat_sat_beam))+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,83,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);%view(10,90); 
%                 hold on; scatter(L1BS.lat_sat_beam(L1BS.beam_ref),(-L1BS.wd_corr(L1BS.beam_ref)+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,80,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','red');
%                 set(gca,'XLim',[min(product_coord(1,product_coord(1,:)~=0)) max(product_coord(1,:))]);
% %                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                 figlabels('Latitude [degrees]','Elevation [m]','','Tracking window position' ,12);
% for i_window=0:L1BS.N_windows
%     hold on; plot(L1BS.lat_sat_beam, zeros(1,length(L1BS.lat_sat_beam))+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
% end
% 
% % subplot(3,3,[5 8]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),coherence_plot'); colormap(colormapMASK);freezeColors;
% %                         figlabels('Beams','Elevation [m]','','Coherence Mask',12);%view(10,26);
% %                         set(gca,'YDir','normal');
% % %                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)])
% % %                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
% coherence_plot=L1BS.coherence_mask;
% coherence_plot(isnan(coherence_plot))=0;
% 
% subplot(3,2,[4 6]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),coherence_plot');
%                         colormap(colormapMASK); freezeColors; set(gca,'YDir','normal');
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
%                         figlabels('Beams','Elevation [m]','', ['Coherence Mask # ' num2str(L1BS.surf_counter, '%03d')] ,12);
%                         for i_window=0:L1BS.N_windows
%                             hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%                         end
%                         
% coherence_plot=L1BS.coherence;
% coherence_plot(isnan(coherence_plot))=0;
% 
% subplot(3,2,[3 5]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),coherence_plot');
%                         colormap(colormapSTACKS); freezeColors; set(gca,'YDir','normal');
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
%                         figlabels('Beams','Elevation [m]','',['Coherence Stack # ' num2str(L1BS.surf_counter, '%03d')] ,12);
%                         for i_window=0:L1BS.N_windows
%                             hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%                         end
%                         colorbar;
% 
%                         
% saveas (h,['Coherence_' num2str(L1BS.surf_counter, '%03d') '.png']);
% close(h);

h=figure; 
colormapPHASE = colormap(jet);
colormapSTACKS=colormap(jet(1024*8));
colormapSTACKS(1,:)=[1 1 1];
% colormapPHASE(1024/2,:)=[1 1 1];
% colormapPHASE(1024/2-1,:)=[1 1 1];
% colormapPHASE(1024/2+1,:)=[1 1 1];

subplot(3,3,1:3); 
                scatter(product_coord(1,:),product_coord(4,:),83,'MarkerEdgeColor',[0 0 0]);
                hold on; scatter(L1BS.lat_sat_beam,(-L1BS.wd_corr(1:length(L1BS.lat_sat_beam))+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,83,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);%view(10,90); 
                hold on; scatter(L1BS.lat_sat_beam(L1BS.beam_ref),(-L1BS.wd_corr(L1BS.beam_ref)+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,80,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','red');
                set(gca,'XLim',[min(product_coord(1,product_coord(1,:)~=0)) max(product_coord(1,:))]);
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
                figlabels('Latitude [degrees]','Elevation [m]','','Tracking window position' ,12);
for i_window=0:L1BS.N_windows
    hold on; plot(L1BS.lat_sat_beam, zeros(1,length(L1BS.lat_sat_beam))+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
end

% subplot(3,3,[5 8]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),coherence_plot'); colormap(colormapMASK);freezeColors;
%                         figlabels('Beams','Elevation [m]','','Coherence Mask',12);%view(10,26);
%                         set(gca,'YDir','normal');
% %                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)])
% %                         set(gca,'XLim',[1 N_max_beams_stack_chd]);

% subplot(3,2,[3 5]); k=surf(1:N_max_beams_stack_chd,y_axis_extend,L1BS.beams_rng_cmpr_Coh');view(10,26);
%                         set(k, 'edgecolor','none');colormap(colormapSTACKS); freezeColors;
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
%                         figlabels('Beams','Elevation [m]','',['Stack # ' num2str(L1BS.surf_counter, '%03d') ' after Coherence masking'] ,12);
%                         for i_window=0:L1BS.N_windows
%                             hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%                         end
%                         colorbar;freezeColors;
% 
% subplot(3,2,[4 6]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),Xdist_filt');
%                         colormap(colormapPHASE); freezeColors; set(gca,'YDir','normal');
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         set(gca,'XLim',[1 N_max_beams_stack_chd]);
%                         figlabels('Beams','Elevation [m]','',['X distance # ' num2str(L1BS.surf_counter, '%03d') ' after Coherence masking'] ,12);
%                         for i_window=0:L1BS.N_windows
%                             hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%                         end
%                         colorbar;freezeColors;
%                         
                        
subplot(3,2,[3 5]); plot(y_axis_extend,nanmean(L1BS.beams_rng_cmpr_Coh)); hold all; plot(y_axis_extend,nanmean(L1BS.beams_rng_cmpr));
                        
                        set(gca,'XLim',[min(y_axis_extend) max(y_axis_extend)]);
                        figlabels('Elevation [m]','FFT p.u','',['Multilooked waveform # ' num2str(L1BS.surf_counter, '%03d')] ,12);
                        for i_window=0:L1BS.N_windows
                            hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
                        end
                        colorbar;freezeColors;

subplot(3,2,[4 6]); plot((y_axis_extend),nanmean(Xdist),'k'); hold all; plot((y_axis_extend),nanmean(Xdist_filt),'g');
                        colormap(colormapPHASE); freezeColors; set(gca,'YDir','normal');
                                                set(gca,'XLim',[min(y_axis_extend) max(y_axis_extend)]);

                        figlabels('Beams','Elevation [m]','',['X distance # ' num2str(L1BS.surf_counter, '%03d')] ,12);
                        for i_window=0:L1BS.N_windows
                            hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
                        end
                        colorbar;freezeColors;
saveas (h,['Beams_and_Xtrack_Extended_Filtered_' num2str(L1BS.surf_counter, '%03d') '.png']);
      
close(h);


end