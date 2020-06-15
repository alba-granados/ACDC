function XStack (L1BS)

global c_cst T0_chd wv_length_ku
global zp_fact_range_cnf N_samples N_max_beams_stack_chd 
global product_coord

%% Axis definition
% 1. Xtrack axis 

Number_beams_Xtrack = 64; % to be added in the config
B = 1.1676;% 1.1676 > Distance between antennas.
deg2xtrack = wv_length_ku * L1BS.win_delay_surf *c_cst/2./(2*pi*B)/180*pi; 

steps=360/Number_beams_Xtrack;
Phase_axis  = (-180:steps:(180));
AoA_axis    = Phase_axis * wv_length_ku./(2*pi*B);
% Xdist_axis  = (Phase_axis-steps/2).*deg2xtrack
Xdist_axis= (-180:360/size(L1BS.phase_diff,2):(180-360/size(L1BS.phase_diff,2))).*deg2xtrack;
% 2. Elevation axis
elev_init = L1BS.alt_sat_beam(L1BS.beam_ref2)-L1BS.win_delay_sar_ku_beam(L1BS.beam_ref2)*c_cst/2;
elev_wd   = L1BS.alt_sat_beam(L1BS.beam_ref)-L1BS.win_delay_sar_ku_beam(L1BS.beam_ref)*c_cst/2; 

% This is the elevation in the nadir direction
nadir_elev_axis = (elev_wd - T0_chd*c_cst/2/zp_fact_range_cnf + N_samples/2*T0_chd*c_cst/2):-T0_chd*c_cst/2/zp_fact_range_cnf:elev_wd-N_samples*T0_chd*c_cst/2+ N_samples/2*T0_chd*c_cst/2;
nadir_elev_axis_extend = (elev_init - T0_chd*c_cst/2/zp_fact_range_cnf + N_samples/2*T0_chd*c_cst/2):-T0_chd*c_cst/2/zp_fact_range_cnf:elev_init-N_samples*L1BS.N_windows*T0_chd*c_cst/2+ N_samples/2*T0_chd*c_cst/2;

% Lets use the different AoA/Phase/Xdist

Slant_Range_axis    = sqrt(Xdist_axis.^2 + (L1BS.win_delay_surf*c_cst/2)^2);
Slant_Elev_axis     = Slant_Range_axis - L1BS.win_delay_surf*c_cst/2 + nadir_elev_axis_extend;
Elev_axis = cos(AoA_axis'.*pi./180)*Slant_Elev_axis;

XStack=nan((length(Phase_axis)-1),size(L1BS.phase_diff,2));
%put NaN values to fill value
phase_diff=L1BS.phase_diff;
phase_diff(isnan(phase_diff))= -999;
for i_sample = 1:size(L1BS.phase_diff,2)
    for i_Xangle=1:(length(Phase_axis)-1)
        
        
%         if(~isempty(phase_values))
            beams_found= find((Phase_axis(i_Xangle)<phase_diff(:,i_sample)).*(phase_diff(:,i_sample)<=Phase_axis(i_Xangle+1)));
            if(~isempty(beams_found))
                XStack(i_Xangle,i_sample)= nanmean(L1BS.beams_rng_cmpr(beams_found,i_sample));
                clear beams_found
            end
           
%         end
        
        
    end
    
end

% figure;k=surf(Slant_Elev_axis,1:N_max_beams_stack_chd, L1BS.beams_rng_cmpr); set(k, 'edgecolor','none');figlabels('Slant Elevation [m]','Beams','','Stack Along track',14); view(250,45); colorbar;
% figure;k=surf(Slant_Elev_axis,1:N_max_beams_stack_chd, L1BS.phase_diff); set(k, 'edgecolor','none');figlabels('Slant Elevation [m]','Beams','','Phase diff Along track',14); colormap(jet);view(-90,90); colorbar;
% figure;k=surf(Slant_Elev_axis, Phase_axis(1:end-1), XStack); set(k, 'edgecolor','none');figlabels('Slant Elevation [m]','Phase diff[deg]','','Stack Xtrack',14);view(70,45); colorbar;
% figure; plot(Phase_axis(1:end-1)*deg2xtrack,nanmean(XStack'));
       

h=figure; 
colormapSTACKS=colormap(jet(1024*8));
colormapSTACKS(1,:)=[1 1 1];

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
% coherence_plot=L1BS.coherence_mask;
% coherence_plot(isnan(coherence_plot))=0;
% subplot(3,3,[5 8]); imagesc(1:N_max_beams_stack_chd,flipud(y_axis_extend),coherence_plot'); colormap(colormapMASK);freezeColors;
%                         figlabels('Beams','Elevation [m]','','Coherence Mask',12);%view(10,26);
%                         set(gca,'YDir','normal');
% %                         set(gca,'YLim',[min(y_axis_extend) max(y_axis_extend)])
% %                         set(gca,'XLim',[1 N_max_beams_stack_chd]);

subplot(3,2,[3 5]); k=surf(1:N_max_beams_stack_chd,nadir_elev_axis_extend,L1BS.beams_rng_cmpr');view(10,26);
                        set(k, 'edgecolor','none');colormap(colormapSTACKS); freezeColors;
                        set(gca,'YLim',[min(nadir_elev_axis_extend) max(nadir_elev_axis_extend)]);
                        set(gca,'XLim',[1 N_max_beams_stack_chd]);
                        figlabels('Beams','Elevation [m]','', ['Stack # ' num2str(L1BS.surf_counter, '%03d') ' Along track'] ,12);
                        for i_window=0:L1BS.N_windows
                            hold on; plot(1:N_max_beams_stack_chd, zeros(1,N_max_beams_stack_chd)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
                        end
        
subplot(3,2,[4 6]); k=surf(Phase_axis(1:end-1)*deg2xtrack,Slant_Elev_axis, XStack');view(10,26);
                        set(k, 'edgecolor','none');colormap(colormapSTACKS); freezeColors;
                        set(gca,'YLim',[min(Slant_Elev_axis) max(Slant_Elev_axis)]);
                        set(gca,'XLim',[min(Phase_axis(1:end-1)*deg2xtrack) max(Phase_axis(1:end-1)*deg2xtrack)]);
                        figlabels('Xtrack distance[m]','Slant Elevation [m]','',['Stack # ' num2str(L1BS.surf_counter, '%03d') ' Across track'],12);
                        resolxtrack = (max(Phase_axis(1:end-1)*deg2xtrack)-min(Phase_axis(1:end-1)*deg2xtrack))/Number_beams_Xtrack;
                        x_axis = min(Phase_axis(1:end-1)*deg2xtrack):resolxtrack:(max(Phase_axis(1:end-1)*deg2xtrack)-resolxtrack);
                        for i_window=0:L1BS.N_windows
                            hold on; plot(x_axis, zeros(1,Number_beams_Xtrack)+elev_init+ N_samples/2*T0_chd*c_cst/2-i_window*N_samples*T0_chd*c_cst/2,'k.')
                        end
%                         
% subplot(2,2,4); plot(y_axis_extend,L1B.wfm_cor_i2q2_sar_ku_Coh,'g');hold all; plot(y_axis_extend,L1B.wfm_cor_i2q2_sar_ku,'b');
%  figlabels('Elevation [m]','FFT power units','',['L1B waveform # ' num2str(L1BS.surf_counter, '%03d')] ,12);
%  set(gca,'XLim',[min(y_axis_extend) max(y_axis_extend)]);
%  legend('Nominal', 'Coherence Filtered')
saveas (h,['Stack_Alongtrack_vs_Xtrack_' num2str(L1BS.surf_counter, '%03d') '.png']);
close(h);


        
        
        
end