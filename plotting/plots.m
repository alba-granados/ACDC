% plots


%%
%Global variables
% global N_total_bursts_sar_ku_isp N_samples_sar_chd N_ku_pulses_burst_chd

%%

% ----FLAGS----
ALL_PLOTS = 1;
ONLY_B1   = 1;
ONLY_B2   = 1;
ONLY_B3   = 1;
ONLY_AUX1 = 1;
SAVE      = 1;
VISIBLE   = 1;

% TEST = [filename_ISP(1:16),filename_ISP(55:56)];


%%

set_default_plot


%% ORBIT TIME OF EACH SURFACE AND SPACING BETWEEN THEM
if ONLY_B1 == 1 || ALL_PLOTS == 1
    %******************************************************
    %**** B1.1 ********************************************
    figure; plot(time_surf,'.')
    set(gca,'XLim',[1 N_total_surf_loc],'FontSize',36) % --> eixos
    xlab = get(gca,'XLabel'); set(xlab,'String','Surface location','FontSize',36)
    ylab = get(gca,'YLabel'); set(ylab,'String','s','FontSize',36)
    if SAVE == 1 && VISIBLE == 0
        figName = [TEST,'_B11_Orbit_time_for_each_surf_loc'];
        print(gcf, figName,'-dpng')
    elseif SAVE == 1 && VISIBLE == 1
        figName = [TEST,'_B11_Orbit_time_for_each_surf_loc'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.fig'])
    end
    
    %******************************************************
    %**** B1.2 ********************************************
    figure; plot(diff(time_surf),'.')
    set(gca,'XLim',[1 N_total_surf_loc-1],'FontSize',36) % --> eixos
    ylab = get(gca,'YLabel'); set(ylab,'String','s','FontSize',36)
    if SAVE == 1
        figName = [TEST,'_B12_Spacing_between_L1_orbit_times'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.fig'])
    end
end


%% STACK DATA
if ONLY_B2 == 1 || ALL_PLOTS == 1
    surf_index = 25:35;
    for i_surf = 1:length(surf_index)
        %******************************************************
        %**** B2.1 ********************************************
        figure; h=surf(1:N_samples_sar_chd*zp_fact_range_cnf,1:N_beams_stack(surf_index(i_surf)),squeeze(beams_rng_cmpr(surf_index(i_surf),1:N_beams_stack(surf_index(i_surf)),:)));
        set(h, 'edgecolor','none');colormap(colormapATDD);view(-30,30)
        set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf],'YLim',[1 N_beams_stack(surf_index(i_surf))])
        xlab = get(gca,'XLabel'); set(xlab,'String','Samples')
        ylab = get(gca,'YLabel'); set(ylab,'String','Beams')
        zlab = get(gca,'ZLabel'); set(zlab,'String','FFT p.u.')
        if SAVE == 1
            figName = [TEST,'_B21_Stack_',num2str(surf_index(i_surf)),'_with_geometry_corr_3D'];
            print(gcf, figName,'-dpng')
            saveas (gcf,[figName,'.fig'])
            saveas (gcf,[figName,'.png'])
        end
        
        %******************************************************
        %**** B2.2 ********************************************
        figure; imagesc(1:N_samples_sar_chd*zp_fact_range_cnf,1:N_beams_stack(surf_index(i_surf)),squeeze(beams_rng_cmpr(surf_index(i_surf),1:N_beams_stack(surf_index(i_surf)),:)))
        set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf],'YLim',[1 N_beams_stack(surf_index(i_surf))])
        xlab = get(gca,'XLabel'); set(xlab,'String','Samples')
        ylab = get(gca,'YLabel'); set(ylab,'String','Beams')
        zlab = get(gca,'ZLabel'); set(zlab,'String','FFT p.u.')
        colorbar('peer',gca,'FontSize',36)
        if SAVE == 1
            figName = [TEST,'_B22_Stack_',num2str(surf_index(i_surf)),'_with_geometry_corr_XY'];
            print(gcf, figName,'-dpng')
            saveas (gcf,[figName,'.fig'])
            saveas (gcf,[figName,'.png'])
        end
    end
end




%% L1 WAVEFORMS
if ONLY_B3 == 1 || ALL_PLOTS == 1
    surf_index = [50,200];
    %******************************************************
    %**** B3.1 ********************************************
    for i_surf = 1:length(surf_index)
        figure; plot(1:N_samples_sar_chd*zp_fact_range_cnf,wfm_cor_i2q2_sar_ku(surf_index(i_surf),:))
        set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf])
        xlab = get(gca,'XLabel'); set(xlab,'String','Samples')
        ylab = get(gca,'YLabel'); set(ylab,'String','FFT p.u.')
        if SAVE == 1
            figName = [TEST,'_B31_L1_wfm_',num2str(surf_index(i_surf))];
            print(gcf, figName,'-dpng')
            saveas (gcf,[figName,'.fig'])
        end
    end
        
    %******************************************************
    %**** B3.2a ********************************************
    figure; imagesc(1:N_samples_sar_chd*zp_fact_range_cnf,1:N_total_surf_loc,wfm_cor_i2q2_sar_ku)
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf],'YLim',[1 N_total_surf_loc])
    xlab = get(gca,'XLabel'); set(xlab,'String','Samples')
    ylab = get(gca,'YLabel'); set(ylab,'String','Stack')
    colorbar('peer',gca,'FontSize',36)
    if SAVE == 1
        figName = [TEST,'_B32a_L1_wfm_XY'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.fig'])
    end
    
    %******************************************************
    %**** B3.2b ********************************************
%     wd_shift    = zeros(1,N_total_surf_loc);
%     wd_shift = ((win_delay_surf-win_delay_surf(1))-(alt_sat-alt_sat(1))*2/c_cst) * bw_ku_chd;
    max_val_l1 = zeros(1,N_total_surf_loc);
    max_pos_l1 = zeros(1,N_total_surf_loc);
%     wfm_cor_i2q2_wdcorr = zeros(N_total_surf_loc,N_samples_sar_chd*zp_fact_range_cnf);
    for i_surf = 1:N_total_surf_loc
%         wd_shift(i_surf) = ((win_delay_surf(i_surf)-win_delay_surf(1))-(alt_sat(i_surf)-alt_sat(1))*2/c_cst) * bw_ku_chd;
%         wfm_cor_i2q2_wdcorr(i_surf,:) = circshift(wfm_cor_i2q2_sar_ku(i_surf,:),[0,round(wd_shift(i_surf)*zp_fact_range_cnf)]);
        [max_val_l1(i_surf),max_pos_l1(i_surf)] = max(wfm_cor_i2q2_sar_ku_wdcorr(i_surf,:));
    end
    figure; imagesc(1:N_samples_sar_chd*zp_fact_range_cnf,1:N_total_surf_loc,wfm_cor_i2q2_sar_ku_wdcorr)
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf],'YLim',[1 N_total_surf_loc])
    xlab = get(gca,'XLabel'); set(xlab,'String','Samples')
    ylab = get(gca,'YLabel'); set(ylab,'String','Stack')
    colorbar('peer',gca,'FontSize',36)
    if SAVE == 1
        figName = [TEST,'_B32b_L1_wfm_XY_aligned'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.fig'])
    end
    
    %******************************************************
    %**** B3.2a 3D ****************************************
    figure; h=surf(1:N_samples_sar_chd*zp_fact_range_cnf,1:N_total_surf_loc,wfm_cor_i2q2_sar_ku);
    set(h, 'edgecolor','none');colorbar;view(-30,30)
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf],'YLim',[1 N_total_surf_loc])
    xlab = get(gca,'XLabel'); set(xlab,'String','Samples')
    ylab = get(gca,'YLabel'); set(ylab,'String','Stack')
    zlab = get(gca,'ZLabel'); set(zlab,'String','FFT p.u.')
    if SAVE == 1
        figName = [TEST,'_B32a_L1_wfm_3D'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.fig'])
        saveas (gcf,[figName,'.png'])
    end
    
    %******************************************************
    %**** B3.2b 3D ****************************************
    figure; h=surf(1:N_samples_sar_chd*zp_fact_range_cnf,1:N_total_surf_loc,wfm_cor_i2q2_sar_ku_wdcorr);
    set(h, 'edgecolor','none');view(-30,30); colorbar;
    set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf],'YLim',[1 N_total_surf_loc])
    xlab = get(gca,'XLabel'); set(xlab,'String','Samples')
    ylab = get(gca,'YLabel'); set(ylab,'String','Stack')
    zlab = get(gca,'ZLabel'); set(zlab,'String','FFT p.u.')
    if SAVE == 1
        figName = [TEST,'_B32b_L1_wfm_3D_aligned'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.fig'])
    end
end












%% ISP WAVEFORMS and ISP WAVEFORMS CAL-GAIN CORRECTED
if ONLY_AUX1 == 1
% 1. Average waveforms within a burst, containing all the data set, with AGC
% attenuations applied.
    zp = 4;
    wfm_isp_av              = zeros(N_total_bursts_sar_ku_isp,N_samples_sar_chd);
    wfm_cal_gain_av         = zeros(N_total_bursts_sar_ku_isp,N_samples_sar_chd*zp);
    wd_shift                = zeros(1,N_total_bursts_sar_ku_isp);
    wfm_cal_gain_wdcorr     = zeros(N_total_bursts_sar_ku_isp,N_samples_sar_chd*zp);
    
    clear wfmAUX wfmAUX2 max_isp max_cal_gain
    for i_burst = 1:N_total_bursts_sar_ku_isp
        % FFT & power
%         wfmAUX = abs(fftshift(fft(squeeze(wfm_iq_sar_ku_isp(i_burst,:,:)).',N_samples_sar_chd*zp),1).').^2/(N_samples_sar_chd);
        wfmAUX2 = abs(fftshift(fft(squeeze(wfm_cal_gain_corrected(i_burst,:,:)).',N_samples_sar_chd*zp),1).').^2/(N_samples_sar_chd);
        
        % Average
        for i_sample = 1:N_samples_sar_chd*zp
%             wfm_isp_av(i_burst,i_sample) = sum(wfmAUX(:,i_sample))/N_ku_pulses_burst_chd;
            wfm_cal_gain_av(i_burst,i_sample) = sum(wfmAUX2(:,i_sample))/N_ku_pulses_burst_chd;
        end
    end
%     max_isp = max(wfm_isp_av.');
    max_cal_gain = max(wfm_cal_gain_av.');
    
    
    %******************************************************
    % ISP waveforms max value
%     figure; plot((1:N_total_bursts_sar_ku_isp),max_isp)
%     hold on
%     plot((1:N_total_bursts_sar_ku_isp),max_cal_gain,'r')
%     set(gca,'XLim',[1 N_total_bursts_sar_ku_isp])
%     xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
%     ylab = get(gca,'YLabel'); set(ylab,'String','FFT p.u.')
%     if SAVE == 1
%         figName = [TEST,'_A11_max_isp_waveforms_with_ATT_and_CAL_applied'];
%         print(gcf, figName,'-dpng')
%         saveas (gcf,[figName,'.fig'])
%     end
    
    
    %******************************************************
    % ISP waveforms
%     figure; mesh(1:N_samples_sar_chd*zp, 1:N_total_bursts_sar_ku_isp, wfm_isp_av)
%     figure; mesh(1:N_samples_sar_chd*zp, 1:N_total_bursts_sar_ku_isp, wfm_cal_gain_av)
%     az = 90; el = 90; view(az,el);
%     figure; imagesc(1:N_samples_sar_chd*zp, 1:N_total_bursts_sar_ku_isp, wfm_isp_av.')
    figure; imagesc(1:N_total_bursts_sar_ku_isp, 1:N_samples_sar_chd*zp, wfm_cal_gain_av.')
    set(gca,'XLim',[1 N_total_bursts_sar_ku_isp])
    set(gca,'YLim',[1 N_samples_sar_chd*zp])
    xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
    ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
    if SAVE == 1
        figName = [TEST,'_A12_isp_waveforms_with_ATT_and_CAL_applied'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.fig'])
    end
    
    
    %******************************************************
    % ISP waveforms with H0 and L1_WinDelay applied
%     wfm_cal_gain_av_aux = [wfm_cal_gain_av, zeros(N_total_bursts_sar_ku_isp,N_samples_sar_chd*zp)];
    wfm_cal_gain_av_aux = wfm_cal_gain_av;
    % adding H0 to FBR waveforms
    for i_burst = 1:N_total_bursts_sar_ku_isp
%         wd_shift(i_burst) = (win_delay_sar_ku(i_burst)-alt_sar_sat(i_burst)*2/c_cst) * bw_ku_chd;
%         wd_shift(i_burst) = win_delay_sar_ku(i_burst) * bw_ku_chd;
        wd_shift(i_burst) = win_delay_sar_ku(i_burst) * bw_ku_chd + 20;
        wfm_cal_gain_wdcorr(i_burst,:) = circshift(wfm_cal_gain_av_aux(i_burst,:),[0,round(wd_shift(i_burst)*zp)]);
    end
%     wfm_cal_gain_wdcorr = fftshift(wfm_cal_gain_wdcorr,2);
    
%     figure; mesh(1:N_samples_sar_chd*zp*2, 1:N_total_bursts_sar_ku_isp, wfm_cal_gain_wdcorr)
%     az = 90; el = 90; view(az,el);
    figure; imagesc(1:N_total_bursts_sar_ku_isp, 1:N_samples_sar_chd*zp, wfm_cal_gain_wdcorr.')
    set(gca,'XLim',[1 N_total_bursts_sar_ku_isp]) % --> eixos
    set(gca,'YLim',[1 N_samples_sar_chd*zp]) % --> eixos
    xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
    ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
    if SAVE == 1
        figName = [TEST,'_A13_isp_waveforms_with_ATT_CAL_H0_WD_applied'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.fig'])
    end

end


