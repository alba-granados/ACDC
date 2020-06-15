%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VALIDATION L1B/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this algorithm is to cross-check the L1B products:
% isardSAT and ESA for CryoSat-2 data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res,L2]=L1B_validation(filename_L1B_ISD,filename_L1B_ESA,filename_L2_ESA,path_results_comparison,surface_aligned,retrack_flag,RETRACKER)
close all;
set_default_plot;
percent_leading_edge=75;
figures_visible=0;

%removal of outliers in ESA SSH retrievals L2
threshold_std=3.0;

if figures_visible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

%% -------------- Create the folder specific for the data -----------------------------
name_file_1B_ISD=strsplit(filename_L1B_ISD,'\');
name_file_1B_ISD=name_file_1B_ISD(end);
name_file_1B_ESA=strsplit(filename_L1B_ESA,'\');
name_file_1B_ESA=name_file_1B_ESA(end);
aux=strsplit(char(name_file_1B_ISD),'.'); aux=char(aux(1));
clear aux;

%-------------- tEXT FILE COMPARISON --------------------------------------
fid = fopen(strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Evaluation.txt'), 'w');
fprintf(fid,'$---------------- Evaluation -----------------------------------$\n');
fprintf(fid,'ISD input file: '); fprintf(fid,'%s\n',char(name_file_1B_ISD));
fprintf(fid,'ESA input file: '); fprintf(fid,'%s\n',char(name_file_1B_ESA));



%% -------------- Constants definition ------------------------------------
semi_minor_axis_cst = 6356752.3142;
mean_sat_alt_cst = 717242;
sec_in_day_cst = 86400;
c_cst = 299792458;
pi_cst = pi;
flat_coeff_cst = 0.00335281066;
earth_radius_cst = 6378137;
semi_major_axis_cst = 6378137;
fs        =   320 * 1e6;
delta_rb_ESA        =   c_cst*0.5/(fs*2);
delta_rb_ISD        =   c_cst*0.5/(fs*double(ncread(filename_L1B_ISD,'zero_padding_l1b_echo_sar_ku')));

%----- Characterizatin Instrument -----------------------------------------
bw_ku_chd               = 320 * 1e6;
pulse_length_chd        = 44.8 * 1e-6;
TBP=bw_ku_chd*pulse_length_chd; % compression gain of chirp pulse

%----- PTR drift information computation ----------------------------------
PTR_power_drift_slope_sec=0.022/(30*24*60*60); %dB/seconds 0.022 dB/month
date_switch_SIRAL=datevec('October 21, 2010 00:00:00','mmmm dd, yyyy HH:MM:SS');
date_acq_end=datevec(ncreadatt(filename_L1B_ISD,'/','last_meas_time'),'dd-mm-yyyy HH:MM:SS');
date_acq_start=datevec(ncreadatt(filename_L1B_ISD,'/','first_meas_time'),'dd-mm-yyyy HH:MM:SS');
%correction over mean time of acquisition: init end of acquisition as not changing
%that much within a data take
PTR_Power_Drift_RX1=(etime(date_acq_start,date_switch_SIRAL)+etime(date_acq_end,date_switch_SIRAL))/2*PTR_power_drift_slope_sec;
%------- Hamming window parameters ----------------------------------------
c1=0.08;
c2=0.92;
N_b=64;
HAM_wind=c1+c2*(cos(pi*(0:1:(N_b-1))/(N_b-1)-pi/2)).^2;
SCALE_HAM=10*log10((mean(HAM_wind)).^2); %dB according to note Dinardo


%% -------------- Read ESA L1B product ------------------------------------
[~,CS]=Cryo_L1b_read(filename_L1B_ESA);
s=size(CS.SAR.data);
ESA_num_surfaces=s(2)*s(3);
ESA_N_samples=s(1);
%---------------- Window Delay --------------------------------------------
ESA_win_delay=reshape(CS.MEA.win_delay,[1,ESA_num_surfaces]);
%---------------- Geometry variables --------------------------------------
ESA_lat_surf=reshape(CS.GEO.LAT,[1,ESA_num_surfaces]);
ESA_lon_surf=reshape(CS.GEO.LON,[1,ESA_num_surfaces]);
ESA_alt_sat=reshape(CS.GEO.H,[1,ESA_num_surfaces]);
ESA_alt_surf=ESA_alt_sat-ESA_win_delay.*c_cst/2;
ESA_xyz_surf=lla2ecef([ESA_lat_surf.',ESA_lon_surf.',ESA_alt_surf.'],flat_coeff_cst,semi_major_axis_cst);         
%----------------- Attitude Information -----------------------------------
ESA_pitch=reshape(CS.GEO.Antenna_Bench_Pitch,[1,ESA_num_surfaces]);
ESA_roll=reshape(CS.GEO.Antenna_Bench_Roll,[1,ESA_num_surfaces]);
ESA_yaw=reshape(CS.GEO.Antenna_Bench_Yaw,[1,ESA_num_surfaces]);

%----------------- Waveforms ----------------------------------------------
ESA_i2q2_meas=reshape(CS.SAR.data,[ESA_N_samples,ESA_num_surfaces]).';
ESA_scale_power=reshape(CS.SAR.echo_scale_power,[ESA_num_surfaces,1]); %known as Scale power: power of 2
ESA_scale_factor=reshape(CS.SAR.echo_scaling,[ESA_num_surfaces,1]); %known as Echo scale factor
%apply scaling factor to waveforms
ESA_i2q2_meas=ESA_i2q2_meas.*repmat((2.^(ESA_scale_power)).*(1e-9*ESA_scale_factor),1,ESA_N_samples);


%% --------------- Read isardSAT L1B product ------------------------------

%---------------- Time ----------------------------------------------------
L2.ISD_time_surf = ncread(filename_L1B_ISD,'time_l1b_echo_sar_ku').';
L2.ISD_UTC_day   = ncread(filename_L1B_ISD,'UTC_day_l1b_echo_sar_ku').';
L2.ISD_UTC_sec   = ncread(filename_L1B_ISD,'UTC_sec_l1b_echo_sar_ku').';


%---------------- Window Delay --------------------------------------------
ISD_win_delay=double(ncread(filename_L1B_ISD,'range_ku_l1b_echo_sar_ku')).';
if surface_aligned
    ISD_win_delay_surf_aligned=double(ncread(filename_L1B_ISD,'range_ku_surf_aligned_l1b_echo_sar_ku')).';
end
ISD_num_surfaces=length(ISD_win_delay);


%----------------- Geophysical corrections --------------------------------
ISD_dry_tropo       = double(ncread(filename_L1B_ISD,'dry_tropo_correction_l1b_echo_sar_ku').');
ISD_wet_tropo       = double(ncread(filename_L1B_ISD,'wet_tropo_correction_l1b_echo_sar_ku').');
ISD_inv_baro        = double(ncread(filename_L1B_ISD,'inverse_baro_correction_l1b_echo_sar_ku').');
ISD_dac             = double(ncread(filename_L1B_ISD,'Dynamic_atmospheric_correction_l1b_echo_sar_ku').');
% ISD_iono_gim        = double(ncread(filename_L1B_ISD,'GIM_iono_correction_l1b_echo_sar_ku').');
ISD_iono_mod        = double(ncread(filename_L1B_ISD,'model_iono_correction_l1b_echo_sar_ku').');
ISD_ocean_tide      = double(ncread(filename_L1B_ISD,'ocean_equilibrium_tide_l1b_echo_sar_ku').');
ISD_lpe_ocean       = double(ncread(filename_L1B_ISD,'long_period_tide_l1b_echo_sar_ku').');
ISD_ocean_loading   = double(ncread(filename_L1B_ISD,'ocean_loading_tide_l1b_echo_sar_ku').');
ISD_solid_earth     = double(ncread(filename_L1B_ISD,'solid_earth_tide_l1b_echo_sar_ku').');
ISD_geoc_polar      = double(ncread(filename_L1B_ISD,'geocentric_polar_tide_l1b_echo_sar_ku').');

ISD_corr_total_sea_ice    = ISD_dry_tropo + ISD_wet_tropo + ISD_inv_baro + ISD_iono_mod+...
                            ISD_ocean_tide + ISD_lpe_ocean + ISD_ocean_loading +...
                            ISD_solid_earth + ISD_geoc_polar;

ISD_corr_total_land_ice   = ISD_dry_tropo + ISD_wet_tropo + ISD_iono_mod +...
                            ISD_ocean_loading + ISD_solid_earth + ISD_geoc_polar;


% The total correction for ocean surfaces has to be retrieved from the L2
% because the Sea State Bias is missing in the FBR file.
% % % ISD_corr_total_ocean      = ISD_dry_tropo + ISD_wet_tropo + ISD_dac + ISD_iono_mod+...
% % %                             ISD_ocean_tide + ISD_lpe_ocean + ISD_ocean_loading +...
% % %                             ISD_solid_earth + ISD_geoc_polar +...
% % %                             ISD_ssb;

ISD_win_delay_corr = ISD_win_delay;% + ISD_corr_total_land_ice;
if surface_aligned
    ISD_win_delay_surf_aligned_corr = ISD_win_delay_surf_aligned + ISD_corr_total_land_ice;
end

%---------------- Geometry variables --------------------------------------
L2.ISD_lat_surf = double(ncread(filename_L1B_ISD,'lat_l1b_echo_sar_ku')).';
L2.ISD_lon_surf = double(ncread(filename_L1B_ISD,'lon_l1b_echo_sar_ku')).';
L2.ISD_alt_sat  = double(ncread(filename_L1B_ISD,'alt_l1b_echo_sar_ku')).';
ISD_alt_surf    = L2.ISD_alt_sat-ISD_win_delay_corr;
if surface_aligned
    ISD_alt_surf_aligned = L2.ISD_alt_sat-ISD_win_delay_surf_aligned_corr;         
    ISD_xyz_surf_aligned = lla2ecef([L2.ISD_lat_surf.',L2.ISD_lon_surf.',ISD_alt_surf_aligned.'],flat_coeff_cst,semi_major_axis_cst);         
end
ISD_xyz_surf = lla2ecef([L2.ISD_lat_surf.',L2.ISD_lon_surf.',ISD_alt_surf.'],flat_coeff_cst,semi_major_axis_cst);  

%write the surfaces on a KML
name_data = strsplit(filename_L1B_ISD,'\');
name_data = char(name_data(end));
%lla2kmlWaveforms_noimage(strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Geolocated_track.kml'),name_data,L2.ISD_lat_surf,L2.ISD_lon_surf,ISD_alt_surf, '.');


%----------------- Attitude Information -----------------------------------
ISD_attitude=double(ncread(filename_L1B_ISD,'satellite_mispointing_l1b_sar_echo_ku')).' / (180/pi_cst * 1e7);
ISD_roll = ISD_attitude(:,2).';
ISD_pitch = ISD_attitude(:,1).';
ISD_yaw = ISD_attitude(:,3).';


%----------------- Waveforms ----------------------------------------------
%not aligned
ISD_i2q2_meas = double(ncread(filename_L1B_ISD,'i2q2_meas_ku_l1b_echo_sar_ku').');
%aligned
if surface_aligned
    ISD_i2q2_meas_wdcorr = double(ncread(filename_L1B_ISD,'i2q2_meas_ku_wdcorr_l1b_echo_sar_ku').');
end
%scaling factor
ISD_waveform_scale = ncread(filename_L1B_ISD,'waveform_scale_factor_l1b_echo_sar_ku');
ISD_N_samples = length(ISD_i2q2_meas(1,:));
%apply scaling factor to waveforms
ISD_i2q2_meas = ISD_i2q2_meas.*repmat(ISD_waveform_scale,1,ISD_N_samples);
if surface_aligned
    ISD_i2q2_meas_wdcorr = ISD_i2q2_meas_wdcorr.*repmat(ISD_waveform_scale,1,ISD_N_samples);
end
%sigma0 scaling factor
ISD_sigma0_scaling = ncread(filename_L1B_ISD,'scale_factor_ku_l1b_echo_sar_ku');



%% ------------------- Read L2 ESA product --------------------------------
[~,CS2] = Cryo_L2_read(filename_L2_ESA);
s = size(CS2.MEA.surf_height_r1_20Hz);
records_db = s(1);
ESA_L2_SSH_r1 = reshape(CS2.MEA.surf_height_r1_20Hz,[1,ESA_num_surfaces]);
ESA_L2_sigma0_r1 = reshape(CS2.MEA.backsc_sig_r1_20Hz,[1,ESA_num_surfaces]);
%undo corrections from ESA
ESA_L2_SSH_r1_nocorr=ESA_L2_SSH_r1+reshape((ones(records_db,1))*(CS2.COR.total_ocean.'),[1,ESA_num_surfaces]);


%% ------------------- FILTERING BY LATITUDE ------------------------------
% LAT_min=max([min(ESA_lat_surf),min(L2.ISD_lat_surf)]);
% LAT_max=min([max(ESA_lat_surf),max(L2.ISD_lat_surf)]);
% ISD_indexes_int=L2.ISD_lat_surf<=LAT_max & L2.ISD_lat_surf>=LAT_min;
% ESA_indexes_int=ESA_lat_surf<=LAT_max & ESA_lat_surf>=LAT_min;
% idx=find(ESA_indexes_int==1);
% ESA_indexes_int(idx+1)=1;
% ESA_valid_surfaces=ESA_alt_sat~=0;
% ISD_indexes_int=L2.ISD_lat_surf<=max(L2.ISD_lat_surf) & L2.ISD_lat_surf>=min(L2.ISD_lat_surf);
% ESA_indexes_int=(ESA_lat_surf<=max(ESA_lat_surf) & ESA_lat_surf>=min(ESA_lat_surf));
% ESA_indexes_int=ESA_indexes_int & ESA_valid_surfaces;

%Forcing the number of surfaces of the ISD and ESA product be the same
%assuming first surface not contemplated in ESA product
ISD_indexes_int=ones(1,ISD_num_surfaces);
ISD_indexes_int(1)=0;
ESA_indexes_int=zeros(1,ESA_num_surfaces);
ESA_indexes_int(1:ISD_num_surfaces-1)=1;
ISD_indexes_int=logical(ISD_indexes_int);
ESA_indexes_int=logical(ESA_indexes_int);

ESA_num_surfaces_filtered=length(find(ESA_indexes_int)==1);
L2.ISD_num_surfaces_filtered=length(find(ISD_indexes_int)==1);

idx_int_ESA=find(ESA_indexes_int);

% % % %-------------------- Surface separations ---------------------------------
% % % ESA_surf_dis=sqrt(sum((diff(ESA_xyz_surf(ESA_indexes_int,:),1,1)).^2,2)); % separation between consecutive surfaces
% % % ISD_surf_dis=sqrt(sum((diff(ISD_xyz_surf(ISD_indexes_int,:),1,1)).^2,2));
% % % if surface_aligned
% % %     ISD_surf_dis_surf_aligned=sqrt(sum((diff(ISD_xyz_surf_aligned(ISD_indexes_int,:),1,1)).^2,2)); 
% % % end



[peak_pow_ESA,idx_max_peak_ESA]=max(ESA_i2q2_meas,[],2);
[peak_pow_ISD,idx_max_peak_ISD]=max(ISD_i2q2_meas,[],2);
L2.idx_int_ISD=find(ISD_indexes_int);
if surface_aligned
    [peak_pow_ISD_wdcorr,idx_max_peak_ISD_wdcorr]=max(ISD_i2q2_meas_wdcorr,[],2);
end

%% ------------------------- Retrack data  --------------------------------
if retrack_flag
    if RETRACKER == 0
    %% -------------------- Simple peak retracker -------------------------------
    %Assuming ocean scenario: 75% of peak can be regarded as a first rough
    %approximation of the leading edge
    %no corrections geophysical corrections will be considered
    %********************* Leading/Epoch estimation ***************************
        for i_wfm=1:ESA_num_surfaces_filtered
            dumm=find(find(ESA_i2q2_meas(idx_int_ESA(i_wfm),:)<=percent_leading_edge/100.0*peak_pow_ESA(idx_int_ESA(i_wfm)))<idx_max_peak_ESA(idx_int_ESA(i_wfm)), 1, 'last' );
            if ~isempty(dumm)
                idx_leading_edge_ESA(i_wfm)=dumm;
            else
                %if there is no leading edge or the waveform has displaced that
                %much from the window to the left select the peak as leading
                %edge
                idx_leading_edge_ESA(i_wfm)=idx_max_peak_ESA(idx_int_ESA(i_wfm));
            end
            ESA_retracking_cor(i_wfm)  =   -(ESA_N_samples/2 - idx_leading_edge_ESA(i_wfm))*delta_rb_ESA;
            ESA_range(i_wfm)           =   ESA_win_delay(idx_int_ESA(i_wfm)).*c_cst/2 + ESA_retracking_cor(i_wfm);
            ESA_SSH(i_wfm)             =   ESA_alt_sat(idx_int_ESA(i_wfm)) - ESA_range(i_wfm);
        end

        for i_wfm=1:L2.ISD_num_surfaces_filtered
            dumm=find(find(ISD_i2q2_meas(L2.idx_int_ISD(i_wfm),:)<=percent_leading_edge/100.0*peak_pow_ISD(L2.idx_int_ISD(i_wfm)))<idx_max_peak_ISD(L2.idx_int_ISD(i_wfm)), 1, 'last' );
            if ~isempty(dumm)
                idx_leading_edge_ISD(i_wfm)=dumm;
            else
                %if there is no leading edge or the waveform has displaced that
                %much from the window to the left select the peak as leading
                %edge
                idx_leading_edge_ISD(i_wfm)=idx_max_peak_ISD(L2.idx_int_ISD(i_wfm));
            end
            ISD_retracking_cor(i_wfm)   =   -(ISD_N_samples/2 - idx_leading_edge_ISD(i_wfm))*delta_rb_ISD;
            L2.ISD_range(i_wfm)         =   ISD_win_delay_corr(L2.idx_int_ISD(i_wfm)) + ISD_retracking_cor(i_wfm);
            L2.ISD_SSH(i_wfm)           =   L2.ISD_alt_sat(L2.idx_int_ISD(i_wfm)) - L2.ISD_range(i_wfm);
            %considering window delay correction/alignment
        %     idx_leading_edge_ISD_wdcorr(i_wfm)=find(find(ISD_i2q2_meas_wdcorr(L2.idx_int_ISD(i_wfm),:)<=percent_leading_edge/100.0*peak_pow_ISD_wdcorr(L2.idx_int_ISD(i_wfm)))<idx_max_peak_ISD_wdcorr(L2.idx_int_ISD(i_wfm)), 1, 'last' );
        %     ISD_retracking_cor_wdcorr(i_wfm)  =   -(ISD_N_samples/2 - idx_leading_edge_ISD_wdcorr(i_wfm))*delta_rb_ISD;
        %     L2.ISD_range_wdcorr(i_wfm)           =   ISD_win_delay_surf_aligned_corr(L2.idx_int_ISD(i_wfm)).*c_cst/2 + ISD_retracking_cor_wdcorr(i_wfm);
        %     L2.ISD_SSH_wdcorr(i_wfm)             =   L2.ISD_alt_sat(L2.idx_int_ISD(i_wfm)) - L2.ISD_range_wdcorr(i_wfm);
            %********************* Sigma0 estimation **********************************
            %sigma0 as peak of waveform
            L2.ISD_sigma0(i_wfm)=10*log10(peak_pow_ISD(L2.idx_int_ISD(i_wfm)))+...
                ISD_sigma0_scaling(L2.idx_int_ISD(i_wfm)); % [dB]

        end
    elseif RETRACKER == 1
    %% -------------------- CoG -------------------------------
        for i_wfm=1:L2.ISD_num_surfaces_filtered
            aux_cog = 0;
            for i_sample = 1:ISD_N_samples
                aux_cog = aux_cog + i_sample*ISD_i2q2_meas(i_sample,i_wfm);
            end
            CoG = aux_cog / sum(ISD_i2q2_meas(:,i_wfm));
        end
    elseif RETRACKER == 2
    %% -------------------- OCoG -------------------------------
    %  -------- From CryoVal-LI technical note -----------------
    %  ---------------------------------------------------------
    % The algorithm used for CryoSat is a variant of that used for the
    % ENVISAT RA2, but adapted and tuned for the CryoSat pulse-limited
    % echo. The OCOG 'amplitude' A is computed as follows, for a
    % configurable sub-window of the echo from bin n1 to n2.
    
        k = 0.3; %pre-set value (as in CryoSat-2 mission over land ice)
        n1 = 35; n2 = ISD_N_samples;
        extra_offset = 3.5;
        
        for i_wfm = 1:L2.ISD_num_surfaces_filtered
            A(i_wfm) = sqrt(sum(ISD_i2q2_meas(n1:n2,i_wfm).^4)/sum(ISD_i2q2_meas(n1:n2,i_wfm).^2));
            
            % Find the range bin i0 in the echo that first crosses the OCOG
            % power threshold kA, where k is a configurable proportion of
            % the amplitude (set to 0.3).
            i0(i_wfm) = find(ISD_i2q2_meas(n1:n2,i_wfm) > k*A(i_wfm),1);
            % If all is well the echo is linearly interpolated to find the
            % exact retracking point t0 (idx_leading_edge_ISD) within the window:
            idx_leading_edge_ISD(i_wfm) = extra_offset + i0(i_wfm)-1+((k*A(i_wfm)-ISD_i2q2_meas(i0(i_wfm)-1,i_wfm))/(ISD_i2q2_meas(i0(i_wfm),i_wfm)-ISD_i2q2_meas(i0(i_wfm)-1,i_wfm)));
        end
    end
end

%% ------------------- COMPARISON -----------------------------------------         
% % % %---------------- Geolocations --------------------------------------------
% % % fprintf(fid,'$---------------- Geolocations --------------------------------------------$\n');
% % % figure; geoshow(ESA_lat_surf(ESA_indexes_int),ESA_lon_surf(ESA_indexes_int),...
% % %     'DisplayType','Point','Marker','o','MarkerEdgeColor','blue'); 
% % % hold on; geoshow(L2.ISD_lat_surf(ISD_indexes_int),L2.ISD_lon_surf(ISD_indexes_int),...
% % %     'DisplayType','Point','Marker','*','MarkerEdgeColor','red'); 
% % % xlabel('Longitude [deg.]'); ylabel('Latitude [deg.]');
% % % title('Comparison Geo-locations L1B: ESA & isardSAT');
% % % legend('L1B-ESA','L1B-ISD');
% % % print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_GEOLOCATIONS_ESA_ISD.png'))
% % % 
% % % if ESA_num_surfaces_filtered == L2.ISD_num_surfaces_filtered
% % %     %compute the geolocation errors
% % %     GEO_loc_errors=sqrt(sum((ESA_xyz_surf(ESA_indexes_int,:)-ISD_xyz_surf(ISD_indexes_int,:)).^2,2));
% % %     if surface_aligned
% % %         GEO_loc_errors_wdcorr=sqrt(sum((ESA_xyz_surf(ESA_indexes_int,:)-ISD_xyz_surf_aligned(ISD_indexes_int,:)).^2,2));
% % %     end    
% % %     figure; plot(1:ESA_num_surfaces_filtered,GEO_loc_errors,'*b'); 
% % %     %hold on; plot(1:L2.ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or'); 
% % %     title('Comparison Geolocation errors: ESA & isardSAT');
% % %     xlabel('Surface Index'); ylabel('Distance [m]');
% % %     legend('ESA-ISD');%,'ESA-ISD (wd. corr.)')
% % %     res.GEO.error_max=max(GEO_loc_errors);
% % %     fprintf(fid,'Maximum geolocation error ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.GEO.error_max);
% % %     res.GEO.mean_error=mean(GEO_loc_errors);
% % %     fprintf(fid,'Mean geolocation error ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.GEO.mean_error);
% % %     res.GEO.std_error=std(GEO_loc_errors);
% % %     fprintf(fid,'Std geolocation error ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.GEO.std_error);
% % %     if surface_aligned
% % %         res.GEO.error_max_wdcorr=max(GEO_loc_errors_wdcorr);
% % %         fprintf(fid,'Maximum geolocation error ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.GEO.error_max_wdcorr);
% % %         res.GEO.mean_error_wdcorr=mean(GEO_loc_errors_wdcorr);
% % %         fprintf(fid,'Mean geolocation error ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.GEO.mean_error_wdcorr);
% % %         res.GEO.std_error_wdcorr=std(GEO_loc_errors_wdcorr);
% % %         fprintf(fid,'Std geolocation error ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.GEO.std_error_wdcorr);
% % %     end
% % % end
% % % 
% % % %----------------- Attitude Information -----------------------------------
% % % fprintf(fid,'$---------------- Attitude Information --------------------------------------------$\n');
% % % figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_roll(ESA_indexes_int),'*b');
% % % hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),ISD_roll(ISD_indexes_int),'or');
% % % xlabel('Latitude [deg.]'); ylabel('Roll [deg]');
% % % title('Comparison Roll L1B: ESA & isardSAT');
% % % legend('L1B-ESA','L1B-ISD');
% % % print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_ROLL_ESA_ISD.png'));
% % % figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_pitch(ESA_indexes_int),'*b');
% % % hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),ISD_pitch(ISD_indexes_int),'or');
% % % xlabel('Latitude [deg.]'); ylabel('Pitch [deg]');
% % % title('Comparison Pitch L1B: ESA & isardSAT');
% % % legend('L1B-ESA','L1B-ISD');
% % % print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_PITCH_ESA_ISD.png'));
% % % figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_yaw(ESA_indexes_int),'*b');
% % % hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),ISD_yaw(ISD_indexes_int),'or');
% % % xlabel('Latitude [deg.]'); ylabel('Yaw [deg]');
% % % title('Comparison Yaw L1B: ESA & isardSAT');
% % % legend('L1B-ESA','L1B-ISD');
% % % print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_YAW_ESA_ISD.png'));
% % % if ESA_num_surfaces_filtered==L2.ISD_num_surfaces_filtered
% % %     res.ATT.roll.max_error=max((ESA_roll(ESA_indexes_int)-ISD_roll(ISD_indexes_int)));
% % %     fprintf(fid,'Maximum error roll ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.roll.max_error);
% % %     res.ATT.roll.mean_error=mean((ESA_roll(ESA_indexes_int)-ISD_roll(ISD_indexes_int)));
% % %     fprintf(fid,'Mean error roll ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.roll.mean_error);  
% % %     res.ATT.roll.std_error=std((ESA_roll(ESA_indexes_int)-ISD_roll(ISD_indexes_int)));
% % %     fprintf(fid,'Std error roll ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.roll.std_error);
% % %     res.ATT.roll.RMSE_error=sqrt(mean((ESA_roll(ESA_indexes_int)-ISD_roll(ISD_indexes_int)).^2));
% % %     fprintf(fid,'RMSE error roll ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.roll.RMSE_error);
% % %     res.ATT.pitch.max_error=max((ESA_pitch(ESA_indexes_int)-ISD_pitch(ISD_indexes_int)));
% % %     fprintf(fid,'Maximum error pitch ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.pitch.max_error);
% % %     res.ATT.pitch.mean_error=mean((ESA_pitch(ESA_indexes_int)-ISD_pitch(ISD_indexes_int)));
% % %     fprintf(fid,'Mean error pitch ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.pitch.mean_error);
% % %     res.ATT.pitch.std_error=std((ESA_pitch(ESA_indexes_int)-ISD_pitch(ISD_indexes_int)));
% % %     fprintf(fid,'Std error pitch ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.pitch.std_error);
% % %     res.ATT.pitch.RMSE_error=sqrt(mean((ESA_pitch(ESA_indexes_int)-ISD_pitch(ISD_indexes_int)).^2));
% % %     fprintf(fid,'RMSE error pitch ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.pitch.RMSE_error);
% % %     res.ATT.yaw.max_error=max((ESA_yaw(ESA_indexes_int)-ISD_yaw(ISD_indexes_int)));
% % %     fprintf(fid,'Maximum error yaw ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.yaw.max_error);
% % %     res.ATT.yaw.mean_error=mean((ESA_yaw(ESA_indexes_int)-ISD_yaw(ISD_indexes_int)));
% % %     fprintf(fid,'Mean error yaw ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.yaw.mean_error); 
% % %     res.ATT.yaw.std_error=std((ESA_yaw(ESA_indexes_int)-ISD_yaw(ISD_indexes_int)));
% % %     fprintf(fid,'Std error yaw ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.yaw.std_error);
% % %     res.ATT.yaw.RMSE_error=sqrt(mean((ESA_yaw(ESA_indexes_int)-ISD_yaw(ISD_indexes_int)).^2));
% % %     fprintf(fid,'RMSE error pitch ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.yaw.RMSE_error);
% % % else
% % %     res.ATT.roll.error_bt_maxs=(max(ESA_roll(ESA_indexes_int))-max(ISD_roll(ISD_indexes_int)));
% % %     fprintf(fid,'Errors between maximums roll ESA-ISD [deg]:' ); fprintf(fid,'%.18g\n',res.ATT.roll.error_bt_maxs);
% % %     res.ATT.roll.error_bt_means=(mean(ESA_roll(ESA_indexes_int))-mean(ISD_roll(ISD_indexes_int)));
% % %     fprintf(fid,'Errors between means roll ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.roll.error_bt_means);
% % %     res.ATT.roll.std_ESA=std(ESA_roll(ESA_indexes_int));
% % %     fprintf(fid,'Std roll ESA [m]: '); fprintf(fid,'%.18g\n',res.ATT.roll.std_ESA);
% % %     res.ATT.roll.std_ISD=std(ISD_roll(ISD_indexes_int));
% % %     fprintf(fid,'Std roll ISD [m]: '); fprintf(fid,'%.18g\n',res.ATT.roll.std_ISD);
% % %     res.ATT.pitch.error_bt_maxs=(max(ESA_pitch(ESA_indexes_int))-max(ISD_pitch(ISD_indexes_int)));
% % %     fprintf(fid,'Errors between maximums pitch ESA-ISD [deg]:' ); fprintf(fid,'%.18g\n',res.ATT.pitch.error_bt_maxs);
% % %     res.ATT.pitch.error_bt_means=(mean(ESA_pitch(ESA_indexes_int))-mean(ISD_pitch(ISD_indexes_int)));
% % %     fprintf(fid,'Errors between means pitch ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.pitch.error_bt_means);
% % %     res.ATT.pitch.std_ESA=std(ESA_pitch(ESA_indexes_int));
% % %     fprintf(fid,'Std pitch ESA [m]: '); fprintf(fid,'%.18g\n',res.ATT.pitch.std_ESA);
% % %     res.ATT.pitch.std_ISD=std(ISD_pitch(ISD_indexes_int));
% % %     fprintf(fid,'Std pitch ISD [m]: '); fprintf(fid,'%.18g\n',res.ATT.pitch.std_ISD);
% % %     res.ATT.yaw.error_bt_maxs=(max(ESA_yaw(ESA_indexes_int))-max(ISD_yaw(ISD_indexes_int)));
% % %     fprintf(fid,'Errors between maximums yaw ESA-ISD [deg]:' ); fprintf(fid,'%.18g\n',res.ATT.yaw.error_bt_maxs);
% % %     res.ATT.yaw.error_bt_means=(mean(ESA_yaw(ESA_indexes_int))-mean(ISD_yaw(ISD_indexes_int)));
% % %     fprintf(fid,'Errors between means yaw ESA-ISD [deg]: '); fprintf(fid,'%.18g\n',res.ATT.yaw.error_bt_means);
% % %     res.ATT.yaw.std_ESA=std(ESA_yaw(ESA_indexes_int));
% % %     fprintf(fid,'Std yaw ESA [m]: '); fprintf(fid,'%.18g\n',res.ATT.yaw.std_ESA);
% % %     res.ATT.yaw.std_ISD=std(ISD_yaw(ISD_indexes_int));
% % %     fprintf(fid,'Std yaw ISD [m]: '); fprintf(fid,'%.18g\n',res.ATT.yaw.std_ISD);
% % % end
% % % 
% % % %--------------- Orbital height -------------------------------------------
% % % fprintf(fid,'$---------------- Orbital height --------------------------------------------$\n');
% % % figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_alt_sat(ESA_indexes_int)/1d3,'*b');
% % % hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),L2.ISD_alt_sat(ISD_indexes_int)/1d3,'or');
% % % xlabel('Latitude [deg.]'); ylabel('Height [Km]');
% % % title('Comparison Orbital Height L1B: ESA & isardSAT');
% % % legend('L1B-ESA','L1B-ISD');
% % % print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_HEIGHT_ESA_ISD.png'));
% % % if ESA_num_surfaces_filtered==L2.ISD_num_surfaces_filtered
% % %     res.ALT.max_error=max((ESA_alt_sat(ESA_indexes_int)-L2.ISD_alt_sat(ISD_indexes_int)));
% % %     fprintf(fid,'Maximum error orbital height ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.ALT.max_error);
% % %     res.ALT.mean_error=mean((ESA_alt_sat(ESA_indexes_int)-L2.ISD_alt_sat(ISD_indexes_int)));
% % %     fprintf(fid,'Mean error orbital height ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.ALT.mean_error);    
% % %     res.ALT.std_error=std((ESA_alt_sat(ESA_indexes_int)-L2.ISD_alt_sat(ISD_indexes_int)));
% % %     fprintf(fid,'Std error orbital height ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.ALT.std_error);
% % %     res.ALT.RMSE_error=sqrt(mean((ESA_alt_sat(ESA_indexes_int)-L2.ISD_alt_sat(ISD_indexes_int)).^2));
% % %     fprintf(fid,'RMSE error orbital height ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.ALT.RMSE_error);
% % % else
% % %     res.ALT.error_bt_maxs=(max(ESA_alt_sat(ESA_indexes_int))-max(L2.ISD_alt_sat(ISD_indexes_int)));
% % %     fprintf(fid,'Errors between maximums orbital height ESA-ISD [m]:' ); fprintf(fid,'%.18g\n',res.ALT.error_bt_maxs);
% % %     res.ALT.error_bt_means=(mean(ESA_alt_sat(ESA_indexes_int))-mean(L2.ISD_alt_sat(ISD_indexes_int)));
% % %     fprintf(fid,'Errors between means orbital height ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.ALT.error_bt_means);
% % %     res.ALT.std_ESA=std(ESA_alt_sat(ESA_indexes_int));
% % %     fprintf(fid,'Std orbital height ESA [m]: '); fprintf(fid,'%.18g\n',res.ALT.std_ESA);
% % %     res.ALT.std_ISD=std(L2.ISD_alt_sat(ISD_indexes_int));
% % %     fprintf(fid,'Std orbital height ISD [m]: '); fprintf(fid,'%.18g\n',res.ALT.std_ISD);
% % % end
% % % 
% % % 
% % % %--------------- Separation between surfaces ------------------------------
% % % fprintf(fid,'$--------------- Separation between surfaces ------------------------------$\n');
% % % figure; plot(1:ESA_num_surfaces_filtered-1,ESA_surf_dis,'*b'); 
% % % hold on; plot(1:L2.ISD_num_surfaces_filtered-1,ISD_surf_dis,'or'); 
% % % if surface_aligned
% % %     hold on; plot(1:L2.ISD_num_surfaces_filtered-1,ISD_surf_dis_surf_aligned,'+g'); 
% % % end
% % % xlabel('Surface Index'); ylabel('Distance [m]');
% % % title('Comparison Surface Separation L1B: ESA & isardSAT');
% % % if surface_aligned
% % %     legend('L1B-ESA','L1B-ISD','L1B-ISD win. del. corr.');
% % % else
% % %     legend('L1B-ESA','L1B-ISD');
% % % end
% % % print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_SURFSEPARATION_ESA_ISD.png'))
% % % 
% % % res.SURF.separation.max_ESA=max(ESA_surf_dis);
% % % fprintf(fid,'Maximum Surface separation [m] ESA: '); fprintf(fid,'%.18g\n',res.SURF.separation.max_ESA);
% % % res.SURF.separation.max_ISD=max(ISD_surf_dis);
% % % fprintf(fid,'Maximum Surface separation [m] ISD: '); fprintf(fid,'%.18g\n',res.SURF.separation.max_ISD);
% % % res.SURF.separation.min_ESA=min(ESA_surf_dis);
% % % fprintf(fid,'Minimum Surface separation [m] ESA: '); fprintf(fid,'%.18g\n',res.SURF.separation.min_ESA);
% % % res.SURF.separation.min_ISD=min(ISD_surf_dis);
% % % fprintf(fid,'Minimum Surface separation [m] ISD: '); fprintf(fid,'%.18g\n',res.SURF.separation.min_ISD);
% % % res.SURF.separation.mean_ESA=mean(ESA_surf_dis);
% % % fprintf(fid,'Mean Surface separation [m] ESA: '); fprintf(fid,'%.18g\n',res.SURF.separation.mean_ESA);
% % % res.SURF.separation.mean_ISD=mean(ISD_surf_dis);
% % % fprintf(fid,'Mean Surface separation [m] ISD: '); fprintf(fid,'%.18g\n',res.SURF.separation.mean_ISD);
% % % res.SURF.separation.std_ESA=std(ESA_surf_dis);
% % % fprintf(fid,'Std Surface separation [m] ESA: '); fprintf(fid,'%.18g\n',res.SURF.separation.std_ESA);
% % % res.SURF.separation.std_ISD=std(ISD_surf_dis);
% % % fprintf(fid,'Std Surface separation [m] ISD: '); fprintf(fid,'%.18g\n',res.SURF.separation.std_ISD);
% % % if surface_aligned
% % %     res.SURF.separation.max_ISD_wdcorr=max(ISD_surf_dis_surf_aligned);
% % %     fprintf(fid,'Maximum Surface separation [m] ISD (wd. corr): '); fprintf(fid,'%.18g\n',res.SURF.separation.max_ISD_wdcorr);
% % %     res.SURF.separation.min_ISD_wdcorr=min(ISD_surf_dis_surf_aligned);
% % %     fprintf(fid,'Minimum Surface separation [m] ISD (wd. corr): '); fprintf(fid,'%.18g\n',res.SURF.separation.min_ISD_wdcorr);
% % %     res.SURF.separation.mean_ISD_wdcorr=mean(ISD_surf_dis_surf_aligned);
% % %     fprintf(fid,'Mean Surface separation [m] ISD (wd. corr.): '); fprintf(fid,'%.18g\n',res.SURF.separation.mean_ISD_wdcorr);
% % %     res.SURF.separation.std_ISD_wdcorr=std(ISD_surf_dis_surf_aligned);
% % %     fprintf(fid,'Std Surface separation [m] ISD (wd. corr.): '); fprintf(fid,'%.18g\n',res.SURF.separation.std_ISD_wdcorr);
% % % end


%---------------- Window Delays --------------------------------------------
fprintf(fid,'$---------------- Window Delays --------------------------------------------$\n');
figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_win_delay(ESA_indexes_int)*c_cst/2/1d3,'-b'); 
hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),ISD_win_delay_corr(ISD_indexes_int)/1d3,'-r'); 
if surface_aligned
    hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),ISD_win_delay_surf_aligned_corr(ISD_indexes_int)*c_cst/2/1d3,'-g'); 
end
xlabel('Latitude [deg.]'); ylabel('Range [Km]');
title('Comparison Range L1B: ESA & isardSAT');
if surface_aligned
    legend('L1B-ESA','L1B-ISD','L1B-ISD win. del. corr.');
else
    legend('L1B-ESA','L1B-ISD');
end
print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_WINDELAYS_ESA_ISD.png'))
if ESA_num_surfaces_filtered==L2.ISD_num_surfaces_filtered
    res.WINDELAY.max_error=max((ESA_win_delay(ESA_indexes_int)-ISD_win_delay_corr(ISD_indexes_int)))*c_cst/2;
    fprintf(fid,'Maximum error window delay ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY.max_error);
    res.WINDELAY.mean_error=mean((ESA_win_delay(ESA_indexes_int)-ISD_win_delay_corr(ISD_indexes_int)))*c_cst/2;
    fprintf(fid,'Mean error window delay ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY.mean_error);
    res.WINDELAY.std_error=std((ESA_win_delay(ESA_indexes_int)-ISD_win_delay_corr(ISD_indexes_int)))*c_cst/2;
    fprintf(fid,'Std error window delay ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY.std_error);
    res.WINDELAY.RMSE_error=sqrt(mean((ESA_win_delay(ESA_indexes_int)-ISD_win_delay_corr(ISD_indexes_int)).^2))*c_cst/2;
    fprintf(fid,'RMSE error window delay ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY.RMSE_error);    
    if surface_aligned
        res.WINDELAY_corr.max_error=max((ESA_win_delay(ESA_indexes_int)-ISD_win_delay_surf_aligned_corr(ISD_indexes_int)))*c_cst/2;
        fprintf(fid,'Maximum error window delay ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY_corr.max_error);
        res.WINDELAY_corr.mean_error=mean((ESA_win_delay(ESA_indexes_int)-ISD_win_delay_surf_aligned_corr(ISD_indexes_int)))*c_cst/2;
        fprintf(fid,'Mean error window delay ESA-ISD (wd. corr.)[m]: '); fprintf(fid,'%.18g\n',res.WINDELAY_corr.mean_error);
        res.WINDELAY_corr.std_error=std((ESA_win_delay(ESA_indexes_int)-ISD_win_delay_surf_aligned_corr(ISD_indexes_int)))*c_cst/2;
        fprintf(fid,'Std error window delay ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY_corr.std_error);
        res.WINDELAY_corr.RMSE_error=sqrt(mean((ESA_win_delay(ESA_indexes_int)-ISD_win_delay_surf_aligned_corr(ISD_indexes_int)).^2))*c_cst/2;
        fprintf(fid,'RMSE error window delay ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY_corr.RMSE_error);
    end  
else
    res.WINDELAY.error_bt_maxs=(max(ESA_win_delay(ESA_indexes_int))-max(ISD_win_delay_corr(ISD_indexes_int)))*c_cst/2;
    fprintf(fid,'Errors between maximums window delay ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY.error_bt_maxs);
    res.WINDELAY.error_bt_means=(mean(ESA_win_delay(ESA_indexes_int))-mean(ISD_win_delay_corr(ISD_indexes_int)))*c_cst/2;
    fprintf(fid,'Errors between means window delay ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY.error_bt_means);
    res.WINDELAY.std_ESA=std(ESA_win_delay(ESA_indexes_int))*c_cst/2;
    fprintf(fid,'Std window delay ESA [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY.std_ESA);
    res.WINDELAY.std_ISD=std(ISD_win_delay_corr(ISD_indexes_int))*c_cst/2;
    fprintf(fid,'Std window delay ISD [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY.std_ISD);
    if surface_aligned
        res.WINDELAY_corr.error_bt_maxs=(max(ESA_win_delay(ESA_indexes_int))-max(ISD_win_delay_surf_aligned_corr(ISD_indexes_int)))*c_cst/2;
        fprintf(fid,'Errors between maximums window delay ESA-ISD (wd. corr) [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY_corr.error_bt_maxs);
        res.WINDELAY_corr.error_bt_means=(mean(ESA_win_delay(ESA_indexes_int))-mean(ISD_win_delay_surf_aligned_corr(ISD_indexes_int)))*c_cst/2;
        fprintf(fid,'Errors between means window delay ESA-ISD (wd. corr) [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY_corr.error_bt_means);
        res.WINDELAY_corr.std_ISD=std(ISD_win_delay_surf_aligned_corr(ISD_indexes_int))*c_cst/2;
        fprintf(fid,'Std window delay ISD (wd. corr) [m]: '); fprintf(fid,'%.18g\n',res.WINDELAY_corr.std_ISD);
    end
end

if retrack_flag
    %---------------- Retracked ranges/SSH ------------------------------------
    % First rough approach: Threshold based % peak
    fprintf(fid,'$---------------- Retracked ranges --------------------------------------------$\n');
    fprintf(fid,'$---------------- Rough estimation --------------------------------------------$\n');
    figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_range/1d3,'-b'); 
    hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),L2.ISD_range/1d3,'or'); 
    %hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),L2.ISD_range_wdcorr/1d3,'+g'); 
    xlabel('Latitude [deg.]'); ylabel('Range [Km]');
    title(strcat('Comparison Retracked Range L1B (',num2str(percent_leading_edge),'% peak retracker): ESA & isardSAT'));
    legend('L1B-ESA','L1B-ISD');%,'L1B-ISD win. del. corr.');
    print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_peakbased_retracked_ranges_ESA_ISD.png'))
    if ESA_num_surfaces_filtered==L2.ISD_num_surfaces_filtered
        res.RETRACKED_RANGE.peak_detector.max_error=max((ESA_range-L2.ISD_range));
        fprintf(fid,'Maximum error retracked range (peak_detector) ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE.peak_detector.max_error);
        res.RETRACKED_RANGE.peak_detector.mean_error=mean((ESA_range-L2.ISD_range));
        fprintf(fid,'Mean error retracked range (peak_detector) ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE.peak_detector.mean_error);
        res.RETRACKED_RANGE.peak_detector.std_error=std((ESA_range-L2.ISD_range));
        fprintf(fid,'Std error retracked range (peak_detector) ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE.peak_detector.std_error);
        res.RETRACKED_RANGE.peak_detector.RMSE_error=sqrt(mean((ESA_range-L2.ISD_range).^2));
        fprintf(fid,'RMSE error retracked range ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE.peak_detector.RMSE_error);
    %     res.RETRACKED_RANGE_corr.peak_detector.max_error=max((ESA_range-L2.ISD_range_wdcorr));
    %     fprintf(fid,'Maximum error retracked range (peak_detector) ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE_corr.peak_detector.max_error);
    %     res.RETRACKED_RANGE_corr.peak_detector.mean_error=mean((ESA_range-L2.ISD_range_wdcorr));
    %     fprintf(fid,'Mean error retracked range (peak_detector) ESA-ISD (wd. corr.)[m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE_corr.peak_detector.mean_error);
    %     res.RETRACKED_RANGE_corr.peak_detector.std_error=std((ESA_range-L2.ISD_range_wdcorr));
    %     fprintf(fid,'Std error retracked range (peak_detector) ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE_corr.peak_detector.std_error);
    %     res.RETRACKED_RANGE_corr.peak_detector.RMSE_error=sqrt(mean((ESA_range-L2.ISD_range_wdcorr).^2));
    %     fprintf(fid,'RMSE error retracked range (peak_detector) ESA-ISD (wd. corr.) [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE_corr.peak_detector.RMSE_error);

    else
        res.RETRACKED_RANGE.peak_detector.error_bt_maxs=abs(max(ESA_range)-max(L2.ISD_range));
        fprintf(fid,'Errors between maximums retracked range (peak_detector) ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE.peak_detector.error_bt_maxs);
        res.RETRACKED_RANGE.peak_detector.error_bt_means=abs(mean(ESA_range)-mean(L2.ISD_range));
        fprintf(fid,'Errors between means retracked range (peak_detector) ESA-ISD [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE.peak_detector.error_bt_means);
        res.RETRACKED_RANGE.peak_detector.std_ESA=std(ESA_range);
        fprintf(fid,'Std retracked range (peak_detector) ESA [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE.peak_detector.std_ESA);
        res.RETRACKED_RANGE.peak_detector.std_ISD=std(L2.ISD_range);
        fprintf(fid,'Std retracked range (peak_detector) ISD [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE.peak_detector.std_ISD);
        res.RETRACKED_RANGE_corr.peak_detector.error_bt_maxs=abs(max(ESA_range)-max(L2.ISD_range));
    %     fprintf(fid,'Errors between maximums retracked range (peak_detector) ESA-ISD (wd. corr) [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE_corr.peak_detector.error_bt_maxs);
    %     res.RETRACKED_RANGE_corr.peak_detector.error_bt_means=abs(mean(ESA_range)-mean(L2.ISD_range_wdcorr));
    %     fprintf(fid,'Errors between means retracked range (peak_detector) ESA-ISD (wd. corr) [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE_corr.peak_detector.error_bt_means);
    %     res.RETRACKED_RANGE_corr.peak_detector.std_ISD=std(L2.ISD_range_wdcorr);
    %     fprintf(fid,'Std retracked range (peak_detector) ISD (wd. corr) [m]: '); fprintf(fid,'%.18g\n',res.RETRACKED_RANGE_corr.peak_detector.std_ISD);
    end
    %******************** SSH *************************************************
    fprintf(fid,'$---------------- SSH --------------------------------------------$\n');
    figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_SSH,'-b'); 
    hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),L2.ISD_SSH,'or'); 
    %hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),L2.ISD_SSH_wdcorr,'+g'); 
    plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1(ESA_indexes_int),'*k');
    plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_nocorr(ESA_indexes_int),'^m')
    xlabel('Latitude [deg.]'); ylabel('SSH [m]');
    title(strcat('Comparison SSH: ESA & isardSAT'));
    legend(strcat('L1B-ESA (Threshold-retracker ',num2str(percent_leading_edge),'%)'),...
        strcat('L1B-ISD (Threshold-retracker ',num2str(percent_leading_edge),'%)'),...
        'L2-ESA','L2-ESA (undoing corrections)');
    print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_peakbased_retracked_SSH_ESA_ISD.png'))
    if ESA_num_surfaces_filtered==L2.ISD_num_surfaces_filtered
        %remove outliers:
        ESA_L2_SSH_r1_filtered=ESA_L2_SSH_r1(ESA_indexes_int);
        ESA_L2_SSH_r1_nocorr_filtered=ESA_L2_SSH_r1_nocorr(ESA_indexes_int);
        idx_outliers_corr=ESA_L2_SSH_r1_filtered>(mean(ESA_L2_SSH_r1_filtered)+threshold_std*std(ESA_L2_SSH_r1_filtered)) | ESA_L2_SSH_r1_filtered<(mean(ESA_L2_SSH_r1_filtered)-threshold_std*std(ESA_L2_SSH_r1_filtered));
        idx_outliers_nocorr=ESA_L2_SSH_r1_nocorr_filtered>(mean(ESA_L2_SSH_r1_nocorr_filtered)+threshold_std*std(ESA_L2_SSH_r1_nocorr_filtered)) | ESA_L2_SSH_r1_nocorr_filtered<(mean(ESA_L2_SSH_r1_nocorr_filtered)-threshold_std*std(ESA_L2_SSH_r1_nocorr_filtered));
        res.SSH.peak_detector.RMSE_error_L1B=sqrt(mean((ESA_SSH-L2.ISD_SSH).^2));
        fprintf(fid,'RMSE error SSH ESA L1B-ISD [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.RMSE_error_L1B);
    %     res.SSH.peak_detector.RMSE_error_L1B_wdcorr=sqrt(mean((ESA_SSH-L2.ISD_SSH_wdcorr).^2));
    %     fprintf(fid,'RMSE error SSH ESA L1B-ISD (wdcorr.) [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.RMSE_error_L1B_wdcorr);
        res.SSH.peak_detector.RMSE_error_L2_corr=sqrt(mean((ESA_L2_SSH_r1_filtered(~idx_outliers_corr)-L2.ISD_SSH(~idx_outliers_corr)).^2));
        fprintf(fid,'RMSE error SSH ESA L2-ISD [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.RMSE_error_L2_corr);
        res.SSH.peak_detector.RMSE_error_L2_nocorr=sqrt(mean((ESA_L2_SSH_r1_nocorr_filtered(~idx_outliers_nocorr)-L2.ISD_SSH(~idx_outliers_nocorr)).^2));
        fprintf(fid,'RMSE error SSH ESA L2(no-corrections)-ISD [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.RMSE_error_L2_nocorr);
    %     res.SSH.peak_detector.RMSE_error_L2_corr_wdcorr=sqrt(mean((ESA_L2_SSH_r1(ESA_indexes_int)-L2.ISD_SSH_wdcorr).^2));
    %     fprintf(fid,'RMSE error SSH ESA L2-ISD (wdcorr.) [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.RMSE_error_L2_corr);
    %     res.SSH.peak_detector.RMSE_error_L2_nocorr_wdcorr=sqrt(mean((ESA_L2_SSH_r1_nocorr(ESA_indexes_int)-L2.ISD_SSH_wdcorr).^2));
    %     fprintf(fid,'RMSE error SSH ESA L2(no-corrections)-ISD (wdcorr.) [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.RMSE_error_L2_nocorr_wdcorr);  
        res.SSH.peak_detector.mean_error_L1B=(mean((ESA_SSH-L2.ISD_SSH)));
        fprintf(fid,'Mean error SSH ESA L1B-ISD [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.mean_error_L1B);
    %     res.SSH.peak_detector.mean_error_L1B_wdcorr=(mean((ESA_SSH-L2.ISD_SSH_wdcorr)));
    %     fprintf(fid,'Mean error SSH ESA L1B-ISD (wdcorr.) [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.mean_error_L1B_wdcorr);
        res.SSH.peak_detector.mean_error_L2_corr=(mean((ESA_L2_SSH_r1_filtered(~idx_outliers_corr)-L2.ISD_SSH(~idx_outliers_corr))));
        fprintf(fid,'Mean error SSH ESA L2-ISD [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.mean_error_L2_corr);
        res.SSH.peak_detector.mean_error_L2_nocorr=(mean((ESA_L2_SSH_r1_nocorr_filtered(~idx_outliers_nocorr)-L2.ISD_SSH(~idx_outliers_nocorr))));
        fprintf(fid,'Mean error SSH ESA L2(no-corrections)-ISD [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.mean_error_L2_nocorr);
    %     res.SSH.peak_detector.mean_error_L2_corr_wdcorr=(mean((ESA_L2_SSH_r1(ESA_indexes_int)-L2.ISD_SSH_wdcorr)));
    %     fprintf(fid,'Mean error SSH ESA L2-ISD (wdcorr.) [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.mean_error_L2_corr);
    %     res.SSH.peak_detector.mean_error_L2_nocorr_wdcorr=(mean((ESA_L2_SSH_r1_nocorr(ESA_indexes_int)-L2.ISD_SSH_wdcorr)));
    %     fprintf(fid,'Mean error SSH ESA L2(no-corrections)-ISD (wdcorr.) [m]: '); fprintf(fid,'%.18g\n',res.SSH.peak_detector.mean_error_L2_nocorr_wdcorr);    
    end

    %******************** sigma0 **********************************************
    fprintf(fid,'$---------------- sigma0 ----------------------------------------$\n');
    figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),'-b'); 
    hold on; plot(L2.ISD_lat_surf(ISD_indexes_int),L2.ISD_sigma0,'-r'); 
    xlabel('Latitude [deg.]'); ylabel('sigma0 [dB]');
    title(strcat('Comparison sigma0: ESA-L2 & isardSAT-peak retracker'));
    legend('L2-ESA','L1B-ISD');
    print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_sigma0_peakretracker_ISD_L2_ESA.png'))

    figure; plot(ESA_L2_sigma0_r1(ESA_indexes_int)-L2.ISD_sigma0,'-r'); 
    xlabel('Record'); ylabel('sigma0 [dB]');
    title(strcat('Difference sigma0: ESA-L2 & isardSAT-peak retracker'));
    print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_sigma0_peakretracker_ISD_L2_ESA_difference.png'))
    if ESA_num_surfaces_filtered==L2.ISD_num_surfaces_filtered
        res.SIGMA0.peak_detector.RMSE_error_L2=sqrt(mean((ESA_L2_sigma0_r1(ESA_indexes_int)-L2.ISD_sigma0).^2));
        fprintf(fid,'RMSE error sigma0 ESA L2 & L1B-ISD [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.peak_detector.RMSE_error_L2);
        res.SIGMA0.peak_detector.mean_error_L2=mean(ESA_L2_sigma0_r1(ESA_indexes_int)-L2.ISD_sigma0);
        fprintf(fid,'Mean error sigma0 ESA L2 & L1B-ISD [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.peak_detector.mean_error_L2);

    end
end

%---------------- Waveforms --------------------------------------------
%****************** ISD ************************************************
%not aligned with window delay corrections
figure; 
[X,Y]=meshgrid(1:1:ISD_N_samples,L2.ISD_lat_surf(ISD_indexes_int));
surf(X,Y,ISD_i2q2_meas(ISD_indexes_int,:)); shading interp; colormap('jet'); 
title('Waveforms record for L1B-ISD')
ylabel('Latitude records [deg]'); xlabel('Samples'); zlabel('Power');
view([45,30]);
print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Waveforms_3D_ISD.png'))
figure; 
imagesc(1:1:ISD_N_samples,L2.ISD_lat_surf(ISD_indexes_int),ISD_i2q2_meas(ISD_indexes_int,:)); 
title('Waveforms record for L1B-ISD')
colormap('jet'); ylabel('Latitude records [deg]'); xlabel('Samples');
colorbar;
print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Waveforms_2D_ISD.png'));


% max_wvfm=max(max(ISD_i2q2_meas(ISD_indexes_int,:)));
% indexes_int=find(ISD_indexes_int==1);
% close all;
% figure;
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% set(fig,'units','pixels','position',[1,1,300,300]);
% %set(0,'defaultFigurePosition',[1,1,300,300])
% for i_surf=1:L2.ISD_num_surfaces_filtered    
%     plot(1:1:ISD_N_samples,ISD_i2q2_meas(indexes_int(i_surf),:)/max_wvfm,'-r');
%     axis([1,ISD_N_samples,0,1.0]);
%     xlabel('Samples'); ylabel('Norm. Power');
%     print('-dpng',strcat(path_results_comparison,'subset\','Waveform_ISD_',num2str(i_surf,'%03d'),'.png'));
% end
% close;
% lla2kmlWaveforms(strcat(path_results_comparison,'Waveforms_ISD_wdcorr.kml'),L2.ISD_lat_surf(ISD_indexes_int),...
%                  L2.ISD_lon_surf(ISD_indexes_int),ISD_alt_surf(ISD_indexes_int), '.', ...
%                  strrep(strcat(path_results_comparison,'subset\'),'\','/'));

% %aligned with window delay corrections
% figure; 
% [X,Y]=meshgrid(1:1:ISD_N_samples,L2.ISD_lat_surf(ISD_indexes_int));
% surf(X,Y,ISD_i2q2_meas_wdcorr(ISD_indexes_int,:)); shading interp; colormap('jet'); 
% title('L1B waveforms record')
% ylabel('Latitude records [deg]'); xlabel('Samples'); zlabel('Power');
% view([45,30]);
% % figure;
% % [X,Y]=meshgrid((L2.ISD_lat_surf(ISD_indexes_int)),fliplr(1:1:ISD_N_samples));
% % surf(X,Y,flip((ISD_i2q2_meas_wdcorr(ISD_indexes_int,:).'),1)/max(max(ISD_i2q2_meas_wdcorr(ISD_indexes_int,:)))); shading interp; colormap('jet'); 
% % title('L1B waveforms record (ROI_1: West Pacific area)')
% % xlabel('Latitude records [deg]'); ylabel('Samples'); zlabel('Power');
% % set(gca,'Ydir','reverse')
% % view([45+180,30]);
% print('-dpng',strcat(path_results_comparison,'Waveforms_wdcorr_3D_ISD.png'));
% figure; 
% imagesc(1:1:ISD_N_samples,L2.ISD_lat_surf(ISD_indexes_int),ISD_i2q2_meas_wdcorr(ISD_indexes_int,:));
% title('Waveforms record for L1B-ISD')
% colormap('jet'); ylabel('Latitude records [deg]'); xlabel('Samples');
% colorbar;
% print('-dpng',strcat(path_results_comparison,'Waveforms_wdcorr_2D_ISD.png'));



%****************** ESA ************************************************
figure; 
[X,Y]=meshgrid(1:1:ESA_N_samples,ESA_lat_surf(ESA_indexes_int));
surf(X,Y,ESA_i2q2_meas(ESA_indexes_int,:)); shading interp; colormap('jet'); 
title('Waveforms record for L1B-ESA')
ylabel('Latitude records [deg]'); xlabel('Samples'); zlabel('Power');
view([45,30]);
print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Waveforms_3D_ESA.png'));
figure; 
imagesc(1:1:ESA_N_samples,ESA_lat_surf(ESA_indexes_int),ESA_i2q2_meas(ESA_indexes_int,:));
title('Waveforms record for L1B-ESA')
colormap('jet'); ylabel('Latitude records [deg]'); xlabel('Samples');
colorbar;
print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Waveforms_2D_ESA.png'));

figure;
plot(ESA_lat_surf(idx_int_ESA),10*log10(peak_pow_ESA(idx_int_ESA)),'-b'); 
hold on; plot(L2.ISD_lat_surf(L2.idx_int_ISD),10*log10(peak_pow_ISD(L2.idx_int_ISD)),'-r'); 
title('Comparison Peak power per record for L1B-ESA')
xlabel('Latitude records [deg]'); ylabel('[dBW]');
legend('L1B-ESA','L1B-ISD');
print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_peakpower_ESA_ISD.png'))

figure;
plot(10*log10(peak_pow_ESA(idx_int_ESA))-10*log10(peak_pow_ISD(L2.idx_int_ISD)),'-r'); 
title('Difference Peak power per record for L1B-ESA')
xlabel('Record'); ylabel('[dBW]');
print('-dpng',strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Comparison_peakpower_ESA_ISD_difference.png'))

fprintf(fid,'$---------------- Power Waveforms --------------------------------------------$\n');
if ESA_num_surfaces_filtered==L2.ISD_num_surfaces_filtered
    res.WVFMS.RMSE_error_peak=sqrt(mean((10*log10(peak_pow_ESA(idx_int_ESA))-10*log10(peak_pow_ISD(L2.idx_int_ISD))).^2));
    fprintf(fid,'RMSE error peak ESA L1B-ISD [dB]: '); fprintf(fid,'%.18g\n',res.WVFMS.RMSE_error_peak);
    res.WVFMS.mean_error_peak=mean(10*log10(peak_pow_ESA(idx_int_ESA))-10*log10(peak_pow_ISD(L2.idx_int_ISD)));
    fprintf(fid,'Mean error peak ESA L1B-ISD [dB]: '); fprintf(fid,'%.18g\n',res.WVFMS.mean_error_peak);
end
res.WVFMS.max_ESA=max(max(ESA_i2q2_meas(ESA_indexes_int,:)));
fprintf(fid,'Maximum value of ESA WVFMs record: '); fprintf(fid,'%.18g\n',res.WVFMS.max_ESA);
res.WVFMS.max_ISD=max(max(ISD_i2q2_meas(ISD_indexes_int,:)));
fprintf(fid,'Maximum value of ISD WVFMs record: '); fprintf(fid,'%.18g\n',res.WVFMS.max_ISD);
res.WVFMS.mean_ESA=mean(mean(ESA_i2q2_meas(ESA_indexes_int,:),1),2);
fprintf(fid,'Mean value of ESA WVFMs record: '); fprintf(fid,'%.18g\n',res.WVFMS.mean_ESA);
res.WVFMS.mean_ISD=mean(mean(ISD_i2q2_meas(ISD_indexes_int,:),1),2);
fprintf(fid,'Mean value of ISD WVFMs record: '); fprintf(fid,'%.18g\n',res.WVFMS.mean_ISD);
if surface_aligned
    res.WVFMS.max_ISD_corr=max(max(ISD_i2q2_meas_wdcorr(ISD_indexes_int,:)));
    fprintf(fid,'Maximum value of ISD (wd. corr.) WVFMs record: '); fprintf(fid,'%.18g\n',res.WVFMS.max_ISD_corr);
    res.WVFMS.mean_ISD_corr=mean(mean(ISD_i2q2_meas_wdcorr(ISD_indexes_int,:),1),2);
    fprintf(fid,'Mean value of ISD (wd. corr.) WVFMs record: '); fprintf(fid,'%.18g\n',res.WVFMS.mean_ISD_corr);
end


fclose(fid);
save(strcat(path_results_comparison,strrep(char(name_file_1B_ISD),'.nc','_'),'Evaluation.mat'),'res');







end