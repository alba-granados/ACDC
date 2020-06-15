function plot_extended_stacks(L1BS,beams_rng_cmpr_extended)

global c_cst T0_chd
global zp_fact_range_cnf N_samples N_bursts_cycle_chd N_ku_pulses_burst_chd
global product_coord

set(0,'defaultFigureVisible','on');
colormapMASK=ones(64,3);
colormapMASK(1,:)=[1 0 0];
colormapMASK(end,:)=[0 1 0];
elev_init=-L1BS.win_delay_sar_ku_beam(L1BS.beam_ref2)*c_cst/2+L1BS.alt_sat_beam(L1BS.beam_ref2);
elev_wd=-L1BS.win_delay_sar_ku_beam(L1BS.beam_ref)*c_cst/2+L1BS.alt_sat_beam(L1BS.beam_ref); 

y_axis = (elev_wd - T0_chd*c_cst/2/zp_fact_range_cnf + N_samples/2*T0_chd*c_cst/2):-T0_chd*c_cst/2/zp_fact_range_cnf:elev_wd-N_samples*T0_chd*c_cst/2+ N_samples/2*T0_chd*c_cst/2;
y_axis_extend = (elev_init - T0_chd*c_cst/2/zp_fact_range_cnf + N_samples/2*T0_chd*c_cst/2):-T0_chd*c_cst/2/zp_fact_range_cnf:elev_init-N_samples*L1BS.N_windows*T0_chd*c_cst/2+ N_samples/2*T0_chd*c_cst/2;

plot_option1=0;
plot_option2=1;
colormapSTACKS=colormap(jet(1024*8));
colormapSTACKS(1,:)=[1 1 1];
if(plot_option1)
h=figure; 



        subplot(3,4,1:3); 
                        
                        scatter(product_coord(1,:),product_coord(4,:),83,'MarkerEdgeColor',[0 0 0]);
                        hold on; scatter(L1BS.lat_sat_beam,(-L1BS.wd_corr(1:length(L1BS.lat_sat_beam))+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,83,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);%view(10,90); 
                        hold on; scatter(L1BS.lat_sat_beam(L1BS.beam_ref),(-L1BS.wd_corr(L1BS.beam_ref)+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,80,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','red');
                        set(gca,'XLim',[min(product_coord(1,product_coord(1,:)~=0)) max(product_coord(1,:))]);
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
                        figlabels('Latitude [degrees]','Elevation [m]','','Tracking window position' ,12);
        for i_window=0:L1BS.N_windows
            hold on; plot(L1BS.lat_sat_beam, zeros(1,length(L1BS.lat_sat_beam))+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
        end
%         subplot(3,4,4); plot_dems(L1BS.lat_sat_beam,L1BS.lat_sat_beam); hold on; 
%                         
%                         ground_track(:,1).Geometry = 'Point';
%                         ground_track(:,1).Lon = L1BS.lat_sat_beam;
%                         ground_track(:,1).Lat = L1BS.lat_sat_beam;
%                         ground_track(:,1).Alt = (-L1BS.wd_corr(1:length(L1BS.lat_sat_beam))+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init;
%                         scatter(L1BS.burst_index,(-L1BS.wd_corr(1:length(L1BS.burst_index))+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,83,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);%view(10,90); 
%                         hold on; scatter(L1BS.burst_index(L1BS.beam_ref),(-L1BS.wd_corr(L1BS.beam_ref)+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,80,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','red');
% %                         set(gca,'XLim',[L1BS.burst_index(1) L1BS.burst_index(1)+N_bursts_cycle_chd*N_ku_pulses_burst_chd-1]);
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
%                         figlabels(' index','Elevation [m]','','Tracking window position' ,12);
%         for i_window=0:L1BS.N_windows
%             hold on; plot(L1BS.burst_index(1):(L1BS.burst_index(1)+N_bursts_cycle_chd*N_ku_pulses_burst_chd-1), zeros(1,N_bursts_cycle_chd*N_ku_pulses_burst_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
%         end

        subplot(3,2,[3 5]); k=surf(1:size(L1BS.beams_rng_cmpr(L1BS.ProcessID_beam==58,:),1),y_axis,L1BS.beams_rng_cmpr(L1BS.ProcessID_beam==58,:)'.*L1BS.good_samples(L1BS.ProcessID_beam==58,:)');view(10,26);
                        set(k, 'edgecolor','none');
                        set(gca,'YLim',[min(y_axis) max(y_axis)])
                        colormap(colormapSTACKS);freezeColors;
                        set(gca,'XLim',[1 N_bursts_cycle_chd*N_ku_pulses_burst_chd]);
                        figlabels('Beams','Elevation [m]','',['Stack #' num2str(L1BS.surf_counter, '%03d') ' before extension - 1 range window'],12);
                        hold on; plot(1:N_bursts_cycle_chd*N_ku_pulses_burst_chd, zeros(1,N_bursts_cycle_chd*N_ku_pulses_burst_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2                           ,'k.')
                        hold on; plot(1:N_bursts_cycle_chd*N_ku_pulses_burst_chd, zeros(1,N_bursts_cycle_chd*N_ku_pulses_burst_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-1*N_samples*T0_chd*c_cst/2,'k.')


        subplot(3,2,[4 6]); k=surf(1:size(L1BS.good_samples(L1BS.ProcessID_beam==58,:),1),y_axis,L1BS.good_samples(L1BS.ProcessID_beam==58,:)'); colormap(colormapMASK);freezeColors;view(10,26);
                        set(k, 'edgecolor','none');figlabels('Beams','Elevation [m]','','Geometry Mask',12);
                        set(gca,'YLim',[min(y_axis) max(y_axis)])
                        set(gca,'XLim',[1 N_bursts_cycle_chd*N_ku_pulses_burst_chd]);  
                        set(gca,'ZLim',[0 100])
saveas (h,['Stack_' num2str(L1BS.surf_counter, '%03d') '.png']);

        subplot(3,2,[4 6]); k=surf(1:size(beams_rng_cmpr_extended(L1BS.ProcessID_beam==58,:),1),y_axis_extend,beams_rng_cmpr_extended(L1BS.ProcessID_beam==58,:)');view(10,26);
                        set(k, 'edgecolor','none');colormap(colormapSTACKS); freezeColors;
                        set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
                        set(gca,'XLim',[1 N_bursts_cycle_chd*N_ku_pulses_burst_chd]);
                        figlabels('Beams','Elevation [m]','',['Stack #' num2str(L1BS.surf_counter, '%03d') '  after extension - ' num2str(L1BS.N_windows) ' range windows'] ,12);
                    for i_window=0:L1BS.N_windows
                        hold on; plot(1:size(beams_rng_cmpr_extended(L1BS.ProcessID_beam==58,:),1), zeros(1,size(beams_rng_cmpr_extended(L1BS.ProcessID_beam==58,:),1))+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
                    end
        
saveas (h,['Stack_Extended_' num2str(L1BS.surf_counter, '%03d') '.png']);


elseif(plot_option2)
    
    
    h=figure;
    
                   
                        scatter(product_coord(1,:),product_coord(4,:),83,'MarkerEdgeColor',[0 0 0]);
                        hold on; scatter(L1BS.lat_sat_beam,(-L1BS.wd_corr(1:length(L1BS.lat_sat_beam))+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,83,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);%view(10,90); 
                        hold on; scatter(L1BS.lat_sat_beam(L1BS.beam_ref),(-L1BS.wd_corr(L1BS.beam_ref)+L1BS.wd_corr(L1BS.beam_ref2))*T0_chd*c_cst/2+elev_init,80,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','red');
                        set(gca,'XLim',[min(product_coord(1,product_coord(1,:)~=0)) max(product_coord(1,:))]);
%                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
                        figlabels('Latitude [degrees]','Elevation [m]','','Tracking window position' ,12);
        for i_window=0:L1BS.N_windows
            hold on; plot(L1BS.lat_sat_beam, zeros(1,length(L1BS.lat_sat_beam))+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
        end
        hold off;
            saveas (h,['window delay_' num2str(L1BS.surf_counter, '%03d') '.png']);

    k=surf(1:size(L1BS.beams_rng_cmpr(L1BS.ProcessID_beam==58,:),1),y_axis,L1BS.beams_rng_cmpr(L1BS.ProcessID_beam==58,:)'.*L1BS.good_samples(L1BS.ProcessID_beam==58,:)');view(10,26);
    set(k, 'edgecolor','none');
    set(gca,'YLim',[min(y_axis) max(y_axis)])
    colormap(colormapSTACKS);freezeColors;
    set(gca,'XLim',[1 N_bursts_cycle_chd*N_ku_pulses_burst_chd]);
    figlabels('Beams','Elevation [m]','',['Stack #' num2str(L1BS.surf_counter, '%03d') ' before extension - 1 range window'],12);
    hold on; plot(1:N_bursts_cycle_chd*N_ku_pulses_burst_chd, zeros(1,N_bursts_cycle_chd*N_ku_pulses_burst_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2                           ,'k.')
    hold on; plot(1:N_bursts_cycle_chd*N_ku_pulses_burst_chd, zeros(1,N_bursts_cycle_chd*N_ku_pulses_burst_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-1*N_samples*T0_chd*c_cst/2,'k.')
    hold off;
    saveas (h,['Stack_' num2str(L1BS.surf_counter, '%03d') '.png']);

    k=surf(1:size(beams_rng_cmpr_extended(L1BS.ProcessID_beam==58,:),1),y_axis_extend,beams_rng_cmpr_extended(L1BS.ProcessID_beam==58,:)');view(10,26);
    set(k, 'edgecolor','none');colormap(colormapSTACKS); freezeColors;
    set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)]);
    set(gca,'XLim',[1 N_bursts_cycle_chd*N_ku_pulses_burst_chd]);
    figlabels('Beams','Elevation [m]','',['Stack #' num2str(L1BS.surf_counter, '%03d') '  after extension - ' num2str(L1BS.N_windows) ' range windows'] ,12);
    for i_window=0:L1BS.N_windows
        hold on; plot(1:size(beams_rng_cmpr_extended(L1BS.ProcessID_beam==58,:),1), zeros(1,size(beams_rng_cmpr_extended(L1BS.ProcessID_beam==58,:),1))+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
    end
    
    saveas (h,['Stack_Extended_' num2str(L1BS.surf_counter, '%03d') '.png']);
    
end  
close(h);
end