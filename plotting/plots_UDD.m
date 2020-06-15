%plots UDD

clear
clc

c_cst = 299792458;
SAVE = 1;
VISIBLE = 1;
set_default_plot;
zp_fact_range_cnf=2;


% PATH = 'C:\Users\roger.ISARDSAT\Documents\Feina\Jason-CS\Inputs_Outputs_and_Plots\input_output_data\GPP_Ph4\';
PATH = 'J:\Jason-CS\Inputs_Outputs_and_Plots\input_output_data\GPP_Ph4\';
TEST = 'T02'; %T01, T02 or T03
L1B_LRM = 0;
L1A_HR = 0;
L1BS_HR = 0;
L1B_HR = 0;

switch TEST
    case 'T01'
        DATE = '20190119T064000_20190119T064019';
    case 'T02'
        DATE = '20210929T064000_20210929T064019';
    case 'T03'
        DATE = '20210929T064000_20210929T064019';
end

%%
% % plotting Burst After Azimuth processing
% % execute inside azimuth processing as beams_focused_shifted is not an output
% MaxLevel= 1.1*max(max(abs(fftshift(fft(squeeze(beams_focused_shifted(i_burst,:,:)).',zp_fact_range_cnf * N_samples_sar_chd),1)))/sqrt(zp_fact_range_cnf * N_samples_sar_chd));
% h=figure; k=surf(1:N_beams_sar_ku(i_burst),1/zp_fact_range_cnf:1/zp_fact_range_cnf:N_samples_sar_chd,abs(fftshift(fft(squeeze(beams_focused_shifted(i_burst,:,:)).',zp_fact_range_cnf * N_samples_sar_chd),1))/sqrt(zp_fact_range_cnf * N_samples_sar_chd));
% figlabels('Beams','Samples','','',12);colormap(colormapATDD);colorbar;
% set(k, 'edgecolor','none');view(-20,80);
% set(gca,'YLim',[1 N_samples_sar_chd]);set(gca,'XLim',[1 N_beams_sar_ku(i_burst)]);set(gca,'ZLim',[0 MaxLevel]);
% title(['Burst ' num2str(i_burst, '%02d')]);
% saveas (h, ['./results/' TEST(1:6) '/3D_2_Burst_Focused_',num2str(i_burst, '%04d'),'_in_frequency_domain.png']);
% view(0,90);
% saveas (h, ['./results/' TEST(1:6) '/2D_2_Burst_Focused_',num2str(i_burst, '%04d'),'_in_frequency_domain.png']);
% close (h);


%% plotting L1B LR
if L1B_LRM == 1
    MODE = 'LRM';
    switch TEST
        case 'T01'
            L1B_LRM = readanyNETCDF_V2([PATH,TEST,'\lrm_l1b_product_C1A_v2.nc']);
        case 'T02'
            L1B_LRM = readanyNETCDF_V2([PATH,TEST,'\lrm_l1b_product_C1B_v2.nc']);
        case 'T03'
            L1B_LRM = readanyNETCDF_V2([PATH,TEST,'\lrm_l1b_product_C1C_v2.nc']);
    end

    N_samples = size(L1B_LRM.data.lrm_power_waveform_ku,1);
    N_records = size(L1B_LRM.data.lrm_power_waveform_ku,2);
    wfm.L1B.(MODE).Ku = zeros(N_records,N_samples);
    wfm.L1B.(MODE).C = zeros(N_records,N_samples);
    for i_record = 1:N_records
        wfm.L1B.(MODE).Ku(i_record,:) = double(L1B_LRM.data.lrm_power_waveform_ku(:,i_record))*double(L1B_LRM.data.waveform_scale_factor_ku(i_record));
        wfm.L1B.(MODE).C(i_record,:) = double(L1B_LRM.data.lrm_power_waveform_c(:,i_record))*double(L1B_LRM.data.waveform_scale_factor_c(i_record));
    end
    MaxLevel= 1.1*max(max(wfm.L1B.(MODE).Ku));

    %-- aligning the waveforms before plotting --%
    alt = (double(L1B_LRM.data.com_altitude_ku)/1e-4 + 1.3e6) / c_cst * 395e6;
    range = (double(L1B_LRM.data.altimeter_range_calibrated_ku)/1e-4 + 1.3e6) / c_cst * 395e6; %1-way range at sample Ns/2
    shift = zeros(1,N_records);
    wfm1 = zeros(N_records,N_samples);
    for i_rec = 1:N_records
        shift(i_rec) = (range(i_rec) - range(1)) - (alt(i_rec) - alt(1));
        wfm1(i_rec,:) = circshift(wfm.L1B.(MODE).Ku(i_rec,:),[0,round(shift(i_rec))]);
    end

    h=figure; k=surf(1:N_records,1:N_samples,wfm1');
    figlabels('L1B records','Samples','','',16);colormap(colormapLRM);colorbar;
    set(k, 'edgecolor','none');view(-20,80);
    set(gca,'YLim',[1 N_samples]);set(gca,'XLim',[1 N_records]);set(gca,'ZLim',[0 MaxLevel]);
    % title('L1B waveforms');
    if SAVE == 1
        saveas (h, ['./results/',TEST,'/',MODE,'/3D_0_L1B_aligned_Ku.png']);
        view(0,-90);saveas (h, ['./results/',TEST,'/',MODE,'/2D_0_L1B_aligned_Ku.png']); close (h);
    end

    MaxLevel= 1.1*max(max(wfm.L1B.(MODE).C));
    h=figure; k=surf(1:N_records,1:N_samples,wfm.L1B.(MODE).C');
    figlabels('L1B records','Samples','','',16);colormap(colormapLRM);colorbar;
    set(k, 'edgecolor','none');view(-20,80);
    set(gca,'YLim',[1 N_samples]);set(gca,'XLim',[1 N_records]);set(gca,'ZLim',[0 MaxLevel]);
    if SAVE == 1
        saveas (h, ['./results/',TEST,'/',MODE,'/3D_0_L1B_aligned_C.png']);
        view(0,-90);saveas (h, ['./results/',TEST,'/',MODE '/2D_0_L1B_aligned_C.png']); close (h);
    end

    % save(['./results/',TEST,'/wfm_L1B_LRMK',],'wfm_L1B_LRM_Ku')
    % save(['./results/',TEST,'/wfm_L1B_LRMC',],'wfm_L1B_LRM_C')
    % -- It is all saved at the end of the script --
    
    clear L1B_LRM
end



%% plotting L1A HR
if L1A_HR == 1
    for i_mode = 1:2
        if i_mode == 1
            MODE = 'HRWK';
        elseif i_mode == 2
            MODE = 'HRMK';
        end
        PROD = 'L1A_';
        L1A = readanyNETCDF_V2([PATH,TEST,'\S6_P4_SIM_',MODE,'_',PROD,'_',DATE,'_',TEST,'.nc']);

        N_samples = size(L1A.data.i_samples_ku,1);
        N_pulses = size(L1A.data.i_samples_ku,2);
        N_bursts = size(L1A.data.i_samples_ku,3);
        wfm_L1A = zeros(N_bursts,N_pulses,N_samples);
        wfm_L1A_fft_av = zeros(N_bursts,N_samples*zp_fact_range_cnf);
        wfm_L1A_fft_av_shifted = zeros(N_bursts,N_samples*zp_fact_range_cnf);
        MaxLevelBurst = zeros(1,N_bursts);
        wd_shift = zeros(1,N_bursts);
        for i_burst = 1:N_bursts
            for i_pulse = 1:N_pulses
                wfm_L1A(i_burst,i_pulse,:) = double(L1A.data.i_samples_ku(:,i_pulse,i_burst)) * double(L1A.data.i_scale_factor_ku(i_pulse,i_burst)) +...
                                  1i*double(L1A.data.q_samples_ku(:,i_pulse,i_burst)) * double(L1A.data.q_scale_factor_ku(i_pulse,i_burst));
            end
            wfm_L1A_fft = abs(fftshift(fft(squeeze(wfm_L1A(i_burst,:,:)).',N_samples*zp_fact_range_cnf),1).').^2/(N_samples) * 2; % *2 because there were 3dB missing in the product
            % Average
            for i_sample = 1:N_samples*zp_fact_range_cnf
                wfm_L1A_fft_av(i_burst,i_sample) = sum(wfm_L1A_fft(:,i_sample))/N_pulses;
            end
            MaxLevelBurst(i_burst)= 1.1*max(wfm_L1A_fft_av(i_burst,:));
        end
        MaxLevel= max(MaxLevelBurst);

        clear alt wd wfm2
        
        %-- aligning the waveforms before plotting --%
        alt = (double(L1A.data.com_altitude_ku) * 1e-4 + 1.3e6);
        wd = (double(L1A.data.altimeter_range_calibrated_ku) * 1e-4 + 1.3e6) * 2/c_cst; %window delay at sample Ns/2 
        for i_burst = 1:N_bursts
            wd_shift(i_burst) = (wd(i_burst)-wd(1)-(alt(i_burst)-alt(1))*2/c_cst) * 395e6;
            wfm_L1A_fft_av_shifted(i_burst,:) = circshift(wfm_L1A_fft_av(i_burst,:),[0,round(wd_shift(i_burst)*zp_fact_range_cnf)]);
        end

        elevation_axis = alt(1) - (wd(1) - 128/395e6)*c_cst/2 - (1:512)/zp_fact_range_cnf/395e6*c_cst/2;
        h=figure; k=surf(1:N_bursts,elevation_axis,wfm_L1A_fft_av_shifted.');
        figlabels('Bursts','Elevation [m]','','',16);colormap(colormapATDD);colorbar;
        set(k, 'edgecolor','none');view(-20,80);
        set(gca,'YLim',[min(elevation_axis) max(elevation_axis)]);set(gca,'XLim',[1 N_bursts]);set(gca,'ZLim',[0 MaxLevel]);
        % title('L1A averaged bursts');
        if SAVE == 1
            saveas (h,['./results/',TEST,'/',MODE,'/3D_1_L1A_Bursts_in_frequency_domain.png']);
            view(0,90); saveas (h,['./results/',TEST,'/',MODE,'/2D_1_L1A_Bursts_in_frequency_domain.png']); close (h);
        end

        clear L1A

    end
end



%% L1B-S
if L1BS_HR == 1
    for i_mode = 1:2
        if i_mode == 1
            MODE = 'HRWK';
        elseif i_mode == 2
            MODE = 'HRMK';
        end
        i_surf=200;

        %------- plotting Stack before Mask application -------
        PROD = 'L1BS';
        CASE = '1';
        L1BS_Amb = readanyNETCDF_V2([PATH,TEST,'\S6_P4_SIM_',MODE,'_',PROD,'_',DATE,'_',TEST,CASE,'.nc']);

        N_samples = size(L1BS_Amb.data.look_i_samples_ku,1);
        N_beams = size(L1BS_Amb.data.look_i_samples_ku,2);
        wfm.L1BS.(MODE).Amb = zeros(N_beams,N_samples);
        for i_beam = 1:N_beams
            wfm.L1BS.(MODE).Amb(i_beam,:) = double(L1BS_Amb.data.look_i_samples_ku(:,i_beam,i_surf)) * double(L1BS_Amb.data.i_scale_factor_ku(i_beam,i_surf)) +...
                             1i*double(L1BS_Amb.data.look_q_samples_ku(:,i_beam,i_surf)) * double(L1BS_Amb.data.q_scale_factor_ku(i_beam,i_surf));
        end
        MaxLevel= 1.1*max(max((abs(wfm.L1BS.(MODE).Amb).^2)/2)); % *2 because there are 3dB missing
        h=figure; k=surf(1:N_beams,1:N_samples,(abs(wfm.L1BS.(MODE).Amb').^2)/2); % *2 because there are 3dB missing
        figlabels('Beams','Samples','','',16);colormap(colormapATDD);colorbar;
        set(k, 'edgecolor','none');view(-20,80);
        set(gca,'YLim',[1 N_samples]);set(gca,'XLim',[1 N_beams]);set(gca,'ZLim',[0 MaxLevel]);
        % title(['Stack ' num2str(i_surf, '%02d')]);
        if SAVE == 1
            saveas (h, ['./results/',TEST,'/',MODE,'/3D_2_L1B-S_Amb_Stack_',num2str(i_surf, '%04d'),'_in_frequency_domain.png']);
            view(0,-90);saveas (h,['./results/',TEST,'/',MODE,'/2D_2_L1B-S_Amb_Stack_',num2str(i_surf, '%04d'),'_in_frequency_domain.png']); close (h);
        end

        clear L1BS_Amb



        %------- plotting Stack after Mask application -------
        PROD = 'L1BS';
        CASE = '2';
        L1BS_NoAmb = readanyNETCDF_V2([PATH,TEST,'\S6_P4_SIM_',MODE,'_',PROD,'_',DATE,'_',TEST,CASE,'.nc']);

        N_samples = size(L1BS_NoAmb.data.look_i_samples_ku,1);
        N_beams = size(L1BS_NoAmb.data.look_i_samples_ku,2);
        wfm.L1BS.(MODE).NoAmb = zeros(N_beams,N_samples);
        for i_beam = 1:N_beams
            wfm.L1BS.(MODE).NoAmb(i_beam,:) = double(L1BS_NoAmb.data.look_i_samples_ku(:,i_beam,i_surf)) * double(L1BS_NoAmb.data.i_scale_factor_ku(i_beam,i_surf)) +...
                             1i*double(L1BS_NoAmb.data.look_q_samples_ku(:,i_beam,i_surf)) * double(L1BS_NoAmb.data.q_scale_factor_ku(i_beam,i_surf));
        end
        MaxLevel= 1.1*max(max((abs(wfm.L1BS.(MODE).NoAmb).^2)/2)); % *2 because there are 3dB missing
        h=figure; k=surf(1:N_beams,1:N_samples,(abs(wfm.L1BS.(MODE).NoAmb').^2)/2); % *2 because there are 3dB missing
        figlabels('Beams','Samples','','',16);colormap(colormapATDD);colorbar;
        set(k, 'edgecolor','none');view(-20,80);
        set(gca,'YLim',[1 N_samples]);set(gca,'XLim',[1 N_beams]);set(gca,'ZLim',[0 MaxLevel]);
        % title(['Stack ' num2str(i_surf, '%02d')]);
        if SAVE == 1
            saveas (h, ['./results/' TEST '/' MODE '/3D_3_L1B-S_NoAmb_Stack_',num2str(i_surf, '%04d'),'_in_frequency_domain.png']);
            view(0,-90);saveas (h, ['./results/' TEST '/' MODE '/2D_3_L1B-S_NoAmb_Stack_',num2str(i_surf, '%04d'),'_in_frequency_domain.png']); close (h);
        end

        clear L1BS_NoAmb


    end
end




%% L1B
if L1B_HR == 1
    for i_mode = 1:2
        if i_mode == 1
            MODE = 'HRWK';
        elseif i_mode == 2
            MODE = 'HRMK';
        end
        %% plotting L1B waveforms aligned and with ambiguities
        PROD = 'L1B_';
        CASE = '1';
        L1B_Amb = readanyNETCDF_V2([PATH,TEST,'\S6_P4_SIM_',MODE,'_',PROD,'_',DATE,'_',TEST,CASE,'.nc']);

        % MaxLevel= 1.1*max(max(wfm_cor_i2q2_sar_ku_wdcorr));
        % h=figure; k=surf(1:N_total_surf_loc,1/zp_fact_range_cnf:1/zp_fact_range_cnf:N_samples_sar_chd,wfm_cor_i2q2_sar_ku_wdcorr.');
        % figlabels('L1b records','Samples','','',30);colormap(colormapATDD);colorbar;
        % set(k, 'edgecolor','none');view(-20,80);
        % set(gca,'YLim',[1 N_samples_sar_chd]);set(gca,'XLim',[1 N_total_surf_loc]);set(gca,'ZLim',[0 MaxLevel]);
        % title([TEST(1:2) ' ' TEST(4:6)]);
        % saveas (h, ['./results/' TEST(1:6) '/3D_6_L1b_waveforms_aligned.png']);
        % view(0,90);
        % saveas (h, ['./results/' TEST(1:6) '/2D_6_L1b_waveforms_aligned.png']);close (h);
        N_samples = size(L1B_Amb.data.sar_power_waveform_ku,1);
        N_surf = size(L1B_Amb.data.sar_power_waveform_ku,2);

        wfm.L1B.(MODE).Amb = zeros(N_surf,N_samples);
        for i_surf = 1:N_surf
            wfm.L1B.(MODE).Amb(i_surf,:) = double(L1B_Amb.data.sar_power_waveform_ku(:,i_surf))*double(L1B_Amb.data.waveform_scale_factor_ku(i_surf));
        end
        MaxLevel= 1.1*max(max(wfm.L1B.(MODE).Amb));
        
        alt = (double(L1B_Amb.data.com_altitude_ku) * 1e-4 + 1.3e6);
        wd = (double(L1B_Amb.data.altimeter_range_calibrated_ku) * 1e-4 + 1.3e6) * 2/c_cst; %window delay at sample Ns/2
        wd_shift = zeros(1,N_records);
        wfm2 = zeros(N_records,N_samples*zp_fact_range_cnf);
        for i_surf = 1:N_surf
            wd_shift(i_surf) = ((wd(i_surf)-wd(1))-(alt(i_surf)-alt(1))*2/c_cst) * 395e6;
            wfm2(i_surf,:) = circshift(wfm.L1B.(MODE).Amb(i_surf,:),[0,round(wd_shift(i_surf)*zp_fact_range_cnf)]);
        end
        elevation_axis = alt(1) - (wd(1) - 128/395e6)*c_cst/2 - (1:512)/zp_fact_range_cnf/395e6*c_cst/2;
        h=figure; k=surf(1:N_surf,elevation_axis,wfm2.');
        figlabels('L1B records','Elevation [m]','','',16);colormap(colormapATDD);colorbar;
        set(k, 'edgecolor','none');view(-20,80);
        set(gca,'YLim',[min(elevation_axis) max(elevation_axis)]);set(gca,'XLim',[1 N_surf]);set(gca,'ZLim',[0 MaxLevel]);
        % title('L1B waveforms');
        if SAVE == 1
            saveas (h, ['./results/' TEST '/' MODE '/3D_4_L1B_aligned_Amb.png']);
            view(0,90); saveas (h, ['./results/' TEST '/' MODE '/2D_4_L1B_aligned_Amb.png']); close (h);
        end
        
        % save(['./results/',TEST,'/wfm_L1B_Amb_',MODE],'wfm')
        % -- It is all saved at the end of the script --
        
        clear L1B_Amb




        %% plotting L1B waveforms aligned and without ambiguities
        PROD = 'L1B_';
        CASE = '2';
        L1B_NoAmb = readanyNETCDF_V2([PATH,TEST,'\S6_P4_SIM_',MODE,'_',PROD,'_',DATE,'_',TEST,CASE,'.nc']);
        
        % MaxLevel= 1.1*max(max(wfm_cor_i2q2_sar_ku_wdcorr));
        % h=figure; k=surf(1:N_total_surf_loc,1/zp_fact_range_cnf:1/zp_fact_range_cnf:N_samples_sar_chd,wfm_cor_i2q2_sar_ku_wdcorr.');
        % figlabels('L1b records','Samples','','',30);colormap(colormapATDD);colorbar;
        % set(k, 'edgecolor','none');view(-20,80);
        % set(gca,'YLim',[1 N_samples_sar_chd]);set(gca,'XLim',[1 N_total_surf_loc]);set(gca,'ZLim',[0 MaxLevel]);
        % title([TEST(1:2) ' ' TEST(4:6)]);
        % saveas (h, ['./results/' TEST(1:6) '/3D_6_L1b_waveforms_aligned.png']);
        % view(0,90);
        % saveas (h, ['./results/' TEST(1:6) '/2D_6_L1b_waveforms_aligned.png']);close (h);
        N_samples = size(L1B_NoAmb.data.sar_power_waveform_ku,1);
        N_surf = size(L1B_NoAmb.data.sar_power_waveform_ku,2);
        
        wfm.L1B.(MODE).NoAmb = zeros(N_surf,N_samples);
        for i_surf = 1:N_surf
            wfm.L1B.(MODE).NoAmb(i_surf,:) = double(L1B_NoAmb.data.sar_power_waveform_ku(:,i_surf))*double(L1B_NoAmb.data.waveform_scale_factor_ku(i_surf));
        end
        MaxLevel= 1.1*max(max(wfm.L1B.(MODE).NoAmb));
        
        alt = (double(L1B_NoAmb.data.com_altitude_ku) * 1e-4 + 1.3e6);
        wd = (double(L1B_NoAmb.data.altimeter_range_calibrated_ku) * 1e-4 + 1.3e6) * 2/c_cst; %window delay at sample Ns/2
        wd_shift = zeros(1,N_records);
        wfm3 = zeros(N_records,N_samples*zp_fact_range_cnf);
        for i_surf = 1:N_surf
            wd_shift(i_surf) = ((wd(i_surf)-wd(1))-(alt(i_surf)-alt(1))*2/c_cst) * 395e6;
            wfm3(i_surf,:) = circshift(wfm.L1B.(MODE).NoAmb(i_surf,:),[0,round(wd_shift(i_surf)*zp_fact_range_cnf)]);
        end
        elevation_axis = alt(1) - (wd(1) - 128/395e6)*c_cst/2 - (1:512)/zp_fact_range_cnf/395e6*c_cst/2;
        h=figure; k=surf(1:N_surf,elevation_axis,wfm3.');
        figlabels('L1B records','Elevation [m]','','',16);colormap(colormapATDD);colorbar;
        set(k, 'edgecolor','none');view(-20,80);
        set(gca,'YLim',[min(elevation_axis) max(elevation_axis)]);set(gca,'XLim',[1 N_surf]);set(gca,'ZLim',[0 MaxLevel]);
        % title('L1B waveforms');
        if SAVE == 1
            saveas (h,['./results/',TEST,'/',MODE,'/3D_5_L1B_aligned_NoAmb.png']);
            view(0,90); saveas (h,['./results/',TEST,'/',MODE,'/2D_5_L1B_aligned_NoAmb.png']); close (h);
        end
        
        
        clear L1B_NoAmb

    end
end

%% Saving L1B waveforms
% % % % % save(['./results/',TEST,'/wfm'],'wfm','TEST','i_surf','N_samples')




% clear



if(0)
    %% plotting a comparison of L1 waveforms
    list = dir(['./results/',TEST,'/*wfm*.mat']);
    for i_file = 1:length(list)
        load(['./results/',TEST,'/',list(i_file).name]);
    end
    h=figure;
    plot(1:N_samples,wfm.L1B.SAR.Amb(i_surf,:),'LineWidth',2)
    hold on
    plot(1:N_samples,wfm.L1B.SAR.NoAmb(i_surf,:),'r','LineWidth',2)
    plot(1:N_samples,wfm.L1B.RMC.Amb(i_surf,:),'g','LineWidth',2)
    plot(1:N_samples,wfm.L1B.RMC.NoAmb(i_surf,:),'k','LineWidth',2)
    legend ('L1B RAW Amb','L1B RAW NoAmb','L1B RMC Amb','L1B RMC NoAmb')
    set(gca,'XLim',[1 N_samples]);
    if SAVE == 1
        saveas(h,['./results/',TEST,'/','/L1B_aligned_comparison.png']);
    end
    close(h);
end