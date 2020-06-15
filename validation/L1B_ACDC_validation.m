%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VALIDATION L2/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this algorithm is to cross-check the L2 products:
% isardSAT and ESA for CryoSat-2 data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res]=L1B_ACDC_validation(filename_L1B_ISR,path_results_comparison,varargin)
global generate_kml
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<2 || nargin>(2+6*2))
    error('Wrong number of input parameters');   
end
p = inputParser;
p.addParamValue('filename_L2_ESA',{''},@(x)ischar(x));
p.addParamValue('figures_visible',0);
p.addParamValue('flag_outliers_removal',0);
p.addParamValue('type_outliers_removal','percentiles');
p.addParamValue('sh_name_nc','ssh');
p.addParamValue('geo_mask',[]);
p.parse(varargin{:});
filename_L2_ESA=char(p.Results.filename_L2_ESA);
figures_visible=p.Results.figures_visible;
flag_outliers_removal=p.Results.flag_outliers_removal;
type_outliers_removal=char(p.Results.type_outliers_removal);
sh_name_nc=p.Results.sh_name_nc;
geo_mask=p.Results.geo_mask;
clear p;


close all;
set_default_plot;
%figures_visible=1;

%removal of outliers in ESA SSH retrievals L2
threshold_std=3.0;

%removal of outliers 
%-------using percentile & IQR (Interquartile Range)
% outliers data<(percentil_low-IQR_times*IQR) | data>(percentil_high+IQR_times*IQR)
IQR_times=1.5; %number of IQR 
outlier_percentil_low=25.0;
outlier_percentil_high=75.0;
%--------using hampel filter
hampel_wind=3;% size of window half size
hampel_sigma=3; %number of std deviations to which a sample differ from local median



%sliding window definitions: smoothing function used as fitting
% number of samples at each side of current position
sliding_window_SSH=10;
sliding_window_SWH=10;
sliding_window_sigma0=10;
sliding_window_COR=10;

%Nbins
nbins=250;

%linestyle definition
linestyle_ESA='or';
linestyle_conventional='*b';
linestyle_ACDC='-r';
linestyle_ESA_smooth='-g';
linestyle_ESA_smooth_nocorr='-y';
linestyle_conventional_smooth='-c';
linestyle_ACDC_smooth='-r';
title_name_conventional='Conventional';
title_name_ACDC='ACDC';


if figures_visible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

%% -------------- Create the folder specific for the data -----------------------------
name_file_L1B_ISR=strsplit(filename_L1B_ISR,'/');
name_file_L1B_ISR=char(name_file_L1B_ISR(end));
name_file_L1B_ISR=name_file_L1B_ISR(17:17+30);
if ~isempty(filename_L2_ESA)
    name_file_L2_ESA=strsplit(filename_L2_ESA,'/');
    name_file_L2_ESA=name_file_L2_ESA(end);
end
output_path=strcat(path_results_comparison,name_file_L1B_ISR,'/');%path_results_comparison;
clear aux;
mkdir(output_path);
%-------------- tEXT FILE COMPARISON --------------------------------------
fid = fopen(strcat(output_path,name_file_L1B_ISR,'_Evaluation.txt'), 'w');
fprintf(fid,'$---------------- Evaluation -----------------------------------$\n');
fprintf(fid,'ISR input file: '); fprintf(fid,'%s\n',char(name_file_L1B_ISR));
if ~isempty(filename_L2_ESA)
    fprintf(fid,'ESA input file: '); fprintf(fid,'%s\n',char(name_file_L2_ESA));
end


%% --------------- Read isardSAT L1B product -------------------------------
%---------------- Geometry variables --------------------------------------
ISR_lat_surf=double(ncread(filename_L1B_ISR,'lat_l1b_echo_sar_ku')).';
ISR_lon_surf=double(ncread(filename_L1B_ISR,'lon_l1b_echo_sar_ku')).';
ISR_num_surfaces=length(ISR_lat_surf);
%---------------- Geophysical parameters ----------------------------------
%------------- Conventional retracker -------------------------------------
ISR_L2_SSH_conventional=double(ncread(filename_L1B_ISR,'ssh_ACDC_PRE_20_ku')).';
ISR_L2_SWH_conventional=double(ncread(filename_L1B_ISR,'swh_ACDC_PRE_20_ku')).';
ISR_L2_sigma0_conventional=double(ncread(filename_L1B_ISR,'sig0_ACDC_PRE_20_ku')).';
ISR_L2_COR_conventional=double(ncread(filename_L1B_ISR,'Pearson_corr_ACDC_PRE_20_ku')).';

%------------- ACDC retracker ---------------------------------------------
ISR_L2_SSH_ACDC=double(ncread(filename_L1B_ISR,'ssh_ACDC_20_ku')).';
ISR_L2_SWH_ACDC=double(ncread(filename_L1B_ISR,'swh_ACDC_20_ku')).';
ISR_L2_sigma0_ACDC=double(ncread(filename_L1B_ISR,'sig0_ACDC_20_ku')).';
ISR_L2_COR_ACDC=double(ncread(filename_L1B_ISR,'Pearson_corr_ACDC_20_ku')).';


%write the surfaces on a KML
if generate_kml
    ISR_alt_surf=ISR_lat_surf;
    ISR_alt_surf(:)=0;
    lla2kmlWaveforms_noimage(strcat(output_path,name_file_L1B_ISR,'_Geolocated_track.kml'),name_file_L1B_ISR,ISR_lat_surf(1:100:end),ISR_lon_surf(1:100:end),ISR_alt_surf(1:100:end), '.');
end



%% ------------------- Read L2 ESA product --------------------------------
if ~isempty(filename_L2_ESA)
    [~,CS2]=Cryo_L2_read(filename_L2_ESA);
    s=size(CS2.MEA.surf_height_r1_20Hz);
    records_db=s(2);
    num_bursts_db=s(1);
    ESA_num_surfaces=records_db*num_bursts_db;
    %-------------- Geometry parameters ---------------------------------------
    ESA_lat_surf=reshape(CS2.MEA.LAT_20Hz,[1,ESA_num_surfaces]);
    ESA_lon_surf=reshape(CS2.MEA.LON_20Hz,[1,ESA_num_surfaces]);
    
    %-------------- Geophysical parameters ------------------------------------
    ESA_L2_SSH_r1=reshape(CS2.MEA.surf_height_r1_20Hz,[1,ESA_num_surfaces]);
    ESA_L2_sigma0_r1=reshape(CS2.MEA.backsc_sig_r1_20Hz,[1,ESA_num_surfaces]);
    %undo corrections from ESA
    ESA_L2_SSH_r1_nocorr=ESA_L2_SSH_r1+reshape((ones(num_bursts_db,1))*(CS2.COR.total_ocean.'),[1,ESA_num_surfaces]);
end

%% ------------------- FILTERING BY LATITUDE ------------------------------
%Forcing the number of surfaces of the ISR and ESA product be the same
%assuming first surface not contemplated in ESA product
if ~isempty(filename_L2_ESA)
    idx_not_lat_lon_zeros=~(ESA_lat_surf==0 & ESA_lon_surf==0);
    indices_lat_lon_zeros=find(idx_not_lat_lon_zeros==0);
    if ISR_num_surfaces <=ESA_num_surfaces
        ISR_indexes_int=ones(1,ISR_num_surfaces);
        ISR_indexes_int(1)=0;
        ESA_indexes_int=zeros(1,ESA_num_surfaces);
        ESA_indexes_int(1:ISR_num_surfaces-1)=1;
        for i_index=1:length(indices_lat_lon_zeros)
            if (indices_lat_lon_zeros(i_index)+1)<=ISR_num_surfaces
                ISR_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
            end
        end
        ISR_indexes_int=logical(ISR_indexes_int);
        ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
    else
        %the number of surfaces ESA limits
        ISR_indexes_int=zeros(1,ISR_num_surfaces);
        ISR_indexes_int(2:ESA_num_surfaces+1)=1;
        ISR_indexes_int(1)=0;
        ESA_indexes_int=ones(1,ESA_num_surfaces);
        for i_index=1:length(indices_lat_lon_zeros)
            if (indices_lat_lon_zeros(i_index)+1)<=ISR_num_surfaces
                ISR_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
            end
        end
        ISR_indexes_int=logical(ISR_indexes_int);
        ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
    end
    
    %---------------- geohraphical filtering ----------------------------------
    %reference is alawys ISR
    % ---- Checking whether the track within geomask if available ---------
    if ~isempty(geo_mask)
        ISR_lon_surf_bis=ISR_lon_surf;
        idx_lt_0= ISR_lon_surf_bis<0;
        %longitudes +- values (-180,180)
        if any(idx_lt_0)
            ISR_lon_surf_bis(idx_lt_0)=ISR_lon_surf_bis(idx_lt_0)+360.0;
        end
        
        clear idx_lt_0;
        
        idx_int_geo=inpolygon(ISR_lon_surf_bis,ISR_lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));
        if ~any(idx_int_geo)
            disp(strcat('Track,',{' '},name_file_L2_ISR,{' '},'outside the limits of the geographical mask'))
            return;
        end
        ESA_indexes_int(ESA_indexes_int==1)=idx_int_geo(ISR_indexes_int==1);
        ISR_indexes_int=ISR_indexes_int & idx_int_geo;
    end
    ESA_num_surfaces_filtered=length(find(ESA_indexes_int)==1);
    ISR_num_surfaces_filtered=length(find(ISR_indexes_int)==1);
    
    idx_int_ESA=find(ESA_indexes_int);
    idx_int_ISR=find(ISR_indexes_int);
    
else
    ISR_indexes_int=ones(1,ISR_num_surfaces);
    ISR_indexes_int=logical(ISR_indexes_int);
    ISR_num_surfaces_filtered=length(find(ISR_indexes_int)==1);
    idx_int_ISR=find(ISR_indexes_int);
end

%% --------------------------- APPLY THE GEO CORR ACCORING TO L2 ESA ------
if ~isempty(filename_L2_ESA)
    if ISR_num_surfaces_filtered==ESA_num_surfaces_filtered
        %reads the geophysical correction for L2 ISR and ESA and applies the missing or removes the non-valid corrections as per L2 ESA
        script_geo_corr_a_ESA_L2
        %only valid for the same number of surfaces (forced)
    end
end



%% --------------------------- OUTLIERS COMPUTATION -----------------------
if ~isempty(filename_L2_ESA)
    idx_outliers_ESA_L2_SSH=zeros(1,ESA_num_surfaces_filtered);
    idx_outliers_ESA_L2_sigma0=zeros(1,ESA_num_surfaces_filtered);
end

idx_outliers_ISR_L2_SSH_conventional      = zeros(1,ISR_num_surfaces_filtered);
idx_outliers_ISR_L2_sigma0_conventional   = zeros(1,ISR_num_surfaces_filtered);
idx_outliers_ISR_L2_SWH_conventional      = zeros(1,ISR_num_surfaces_filtered);
idx_outliers_ISR_L2_COR_conventional      = zeros(1,ISR_num_surfaces_filtered);

idx_outliers_ISR_L2_SSH_ACDC      = zeros(1,ISR_num_surfaces_filtered);
idx_outliers_ISR_L2_sigma0_ACDC   = zeros(1,ISR_num_surfaces_filtered);
idx_outliers_ISR_L2_SWH_ACDC      = zeros(1,ISR_num_surfaces_filtered);
idx_outliers_ISR_L2_COR_ACDC      = zeros(1,ISR_num_surfaces_filtered);
    
if flag_outliers_removal
    if ~isempty(filename_L2_ESA)
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ----------------------------- SSH ---------------------------------------
        switch type_outliers_removal
            case 'percentiles'
                [ESA_L2_SSH_r1(ESA_indexes_int),idx_outliers_ESA_L2_SSH]=outliers_by_percentiles(ESA_L2_SSH_r1(ESA_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                ESA_L2_SSH_r1_nocorr(idx_int_ESA(idx_outliers_ESA_L2_SSH))=NaN;
            case 'hampel'
                %[ESA_L2_SSH_r1(ESA_indexes_int)] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1(ESA_indexes_int),hampel_wind,hampel_sigma);
                [~,idx_outliers_ESA_L2_SSH] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1(ESA_indexes_int),hampel_wind,hampel_sigma);
                ESA_L2_SSH_r1(idx_int_ESA(idx_outliers_ESA_L2_SSH))=NaN;
                %[ESA_L2_SSH_r1_nocorr(ESA_indexes_int),indx] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_nocorr(ESA_indexes_int),hampel_wind,hampel_sigma);
                ESA_L2_SSH_r1_nocorr(idx_int_ESA(idx_outliers_ESA_L2_SSH))=NaN;
        end
        % ----------------------------- SIGMA0 ------------------------------------
        switch type_outliers_removal
            case 'percentiles'
                [ESA_L2_sigma0_r1(ESA_indexes_int),idx_outliers_ESA_L2_sigma0]=outliers_by_percentiles(ESA_L2_sigma0_r1(ESA_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
            case 'hampel'
                %[ESA_L2_sigma0_r1(ESA_indexes_int),indx] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),hampel_wind,hampel_sigma);
                [~,idx_outliers_ESA_L2_sigma0] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),hampel_wind,hampel_sigma);
                ESA_L2_sigma0_r1(idx_int_ESA(idx_outliers_ESA_L2_sigma0))=NaN;
        end        
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ISR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------------- Conventional -----------------------------------
    switch type_outliers_removal
        case 'percentiles'
            % compute the errors w.r.t fitting on the data using a smooth
            % function
            %------------------- SSH ----------------------------------
            [ISR_L2_SSH_conventional(ISR_indexes_int),idx_outliers_ISR_L2_SSH_conventional]=outliers_by_percentiles(ISR_L2_SSH_conventional(ISR_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
            
            %----------------- sigma0 ---------------------------------
            [ISR_L2_sigma0_conventional(ISR_indexes_int),idx_outliers_ISR_L2_sigma0_conventional]=outliers_by_percentiles(ISR_L2_sigma0_conventional(ISR_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
            
            %----------------- SWH ------------------------------------
            [ISR_L2_SWH_conventional(ISR_indexes_int),idx_outliers_ISR_L2_SWH_conventional]=outliers_by_percentiles(ISR_L2_SWH_conventional(ISR_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
            
            % ---------- COR: pearson correlation coefficient ---------
            [ISR_L2_COR_conventional(ISR_indexes_int),idx_outliers_ISR_L2_COR_conventional]=outliers_by_percentiles(ISR_L2_COR_conventional(ISR_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
        case 'hampel'
            %------------------ SSH ---------------------------
            %[ISR_L2_SSH_conventional(ISR_indexes_int),indx] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_SSH_conventional(ISR_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ISR_L2_SSH_conventional] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_SSH_conventional(ISR_indexes_int),hampel_wind,hampel_sigma);
            ISR_L2_SSH_conventional(idx_int_ISR(idx_outliers_ISR_L2_SSH_conventional))=NaN;
            %----------------- sigma0 ---------------------------------
            %[ISR_L2_sigma0_conventional(ISR_indexes_int),indx] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_sigma0_conventional(ISR_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ISR_L2_sigma0_conventional] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_sigma0_conventional(ISR_indexes_int),hampel_wind,hampel_sigma);
            ISR_L2_sigma0_conventional(idx_int_ISR(idx_outliers_ISR_L2_sigma0_conventional))=NaN;
            %----------------- SWH ------------------------------------
            %[ISR_L2_SWH_conventional(ISR_indexes_int),indx] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_SWH_conventional(ISR_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ISR_L2_SWH_conventional] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_SWH_conventional(ISR_indexes_int),hampel_wind,hampel_sigma);
            ISR_L2_SWH_conventional(idx_int_ISR(idx_outliers_ISR_L2_SWH_conventional))=NaN;
            % ---------- COR: pearson correlation coefficient ---------
            %[ISR_L2_COR_conventional(ISR_indexes_int),indx] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_COR_conventional(ISR_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ISR_L2_COR_conventional] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_COR_conventional(ISR_indexes_int),hampel_wind,hampel_sigma);
            ISR_L2_COR_conventional(idx_int_ISR(idx_outliers_ISR_L2_COR_conventional))=NaN;
    end
    
    %------------------------ ACDC ----------------------------------------
    switch type_outliers_removal
        case 'percentiles'
            %------------------- SSH ----------------------------------
            [ISR_L2_SSH_ACDC(ISR_indexes_int),idx_outliers_ISR_L2_SSH_ACDC]=outliers_by_percentiles(ISR_L2_SSH_ACDC(ISR_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
            %------------------- SSH ----------------------------------
            [ISR_L2_SWH_ACDC(ISR_indexes_int),idx_outliers_ISR_L2_SWH_ACDC]=outliers_by_percentiles(ISR_L2_SWH_ACDC(ISR_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
            %----------------- sigma0 ---------------------------------
            [ISR_L2_sigma0_ACDC(ISR_indexes_int),idx_outliers_ISR_L2_sigma0_ACDC]=outliers_by_percentiles(ISR_L2_sigma0_ACDC(ISR_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
            % ---------- COR: pearson correlation coefficient ---------
            [ISR_L2_COR_ACDC(ISR_indexes_int),idx_outliers_ISR_L2_COR_ACDC]=outliers_by_percentiles(ISR_L2_COR_ACDC(ISR_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
        case 'hampel'
            %------------------- SSH ----------------------------------
            %[ISR_L2_SSH_ACDC(ISR_indexes_int),indx] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_SSH_ACDC(ISR_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ISR_L2_SSH_ACDC] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_SSH_ACDC(ISR_indexes_int),hampel_wind,hampel_sigma);
            ISR_L2_SSH_ACDC(idx_int_ISR(idx_outliers_ISR_L2_SSH_ACDC))=NaN;
            %------------------- SWH ----------------------------------
            %[ISR_L2_SSH_ACDC(ISR_indexes_int),indx] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_SSH_ACDC(ISR_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ISR_L2_SWH_ACDC] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_SWH_ACDC(ISR_indexes_int),hampel_wind,hampel_sigma);
            ISR_L2_SWH_ACDC(idx_int_ISR(idx_outliers_ISR_L2_SWH_ACDC))=NaN;
            %----------------- sigma0 ---------------------------------
            %[ISR_L2_sigma0_ACDC(ISR_indexes_int),indx] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_sigma0_ACDC(ISR_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ISR_L2_sigma0_ACDC] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_sigma0_ACDC(ISR_indexes_int),hampel_wind,hampel_sigma);
            ISR_L2_sigma0_ACDC(idx_int_ISR(idx_outliers_ISR_L2_sigma0_ACDC))=NaN;
            % ---------- COR: pearson correlation coefficient ---------
            %[ISR_L2_COR_ACDC(ISR_indexes_int),indx] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_COR_ACDC(ISR_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ISR_L2_COR_ACDC] = hampel(ISR_lat_surf(ISR_indexes_int),ISR_L2_COR_ACDC(ISR_indexes_int),hampel_wind,hampel_sigma);
            ISR_L2_COR_ACDC(idx_int_ISR(idx_outliers_ISR_L2_COR_ACDC))=NaN;
    end
    % compute the errors w.r.t fitting on the data using a smooth
    % function

end

%% ------------------- COMPARISON -----------------------------------------         
%--------------------------------------------------------------------------
%------ Compute the errors around a fitting of the SSH & sigma0 -----------
%--------------------------------------------------------------------------
% A smoothing window is used
%--------------------------------------------------------------------------
if ~isempty(filename_L2_ESA)
    %-----------------------------ESA------------------------------------------
    for i_wfm=1:ESA_num_surfaces_filtered
        % compute the errors w.r.t fitting on the data using a smooth
        % function
        %------------------------------ SSH -----------------------------------
        start               = max(i_wfm-sliding_window_SSH,1);
        finish              = min(i_wfm+sliding_window_SSH-1,ESA_num_surfaces_filtered);
        ESA_L2_SSH_r1_smoothed(i_wfm) = nanmean(ESA_L2_SSH_r1(idx_int_ESA(start:finish)));
        ESA_L2_SSH_r1_nocorr_smoothed(i_wfm) = nanmean(ESA_L2_SSH_r1_nocorr(idx_int_ESA(start:finish)));
        
        %-------------------- sigma0 from L2 ESA product ----------------------
        start               = max(i_wfm-sliding_window_sigma0,1);
        finish              = min(i_wfm+sliding_window_sigma0-1,ESA_num_surfaces_filtered);
        ESA_L2_sigma0_r1_smoothed(i_wfm) = nanmean(ESA_L2_sigma0_r1(idx_int_ESA(start:finish)));
    end
end

%--------------------------------------------------------------------------
%--------------------------- ISR ------------------------------------------
%-------------------------- CONVENTIONAL ----------------------------------
ISR_L2_SSH_conventional_smoothed          = zeros(1,ISR_num_surfaces_filtered);
ISR_L2_sigma0_conventional_smoothed       = zeros(1,ISR_num_surfaces_filtered);
ISR_L2_SWH_conventional_smoothed          = zeros(1,ISR_num_surfaces_filtered);
ISR_L2_COR_conventional_smoothed          = zeros(1,ISR_num_surfaces_filtered);
%------------------------------ ACDC --------------------------------------
ISR_L2_SSH_ACDC_smoothed          = zeros(1,ISR_num_surfaces_filtered);
ISR_L2_SWH_ACDC_smoothed          = zeros(1,ISR_num_surfaces_filtered);
ISR_L2_sigma0_ACDC_smoothed       = zeros(1,ISR_num_surfaces_filtered);
ISR_L2_COR_ACDC_smoothed          = zeros(1,ISR_num_surfaces_filtered);

for i_wfm=1:ISR_num_surfaces_filtered
    start_SSH               = max(i_wfm-sliding_window_SSH,1);
    finish_SSH              = min(i_wfm+sliding_window_SSH-1,ISR_num_surfaces_filtered);
    start_sigma0            = max(i_wfm-sliding_window_sigma0,1);
    finish_sigma0           = min(i_wfm+sliding_window_sigma0-1,ISR_num_surfaces_filtered);
    start_SWH               = max(i_wfm-sliding_window_SWH,1);
    finish_SWH              = min(i_wfm+sliding_window_SWH-1,ISR_num_surfaces_filtered);
    start_COR               = max(i_wfm-sliding_window_COR,1);
    finish_COR              = min(i_wfm+sliding_window_COR-1,ISR_num_surfaces_filtered);
    
    %------------------------ CONVENTIONAL --------------------------------
    % compute the errors w.r.t fitting on the data using a smooth
    % function
    %------------------- SSH ----------------------------------------------
    ISR_L2_SSH_conventional_smoothed(i_wfm)   = nanmean(ISR_L2_SSH_conventional(idx_int_ISR(start_SSH:finish_SSH)));
    %----------------- sigma0 ---------------------------------
    ISR_L2_sigma0_conventional_smoothed(i_wfm) = nanmean(ISR_L2_sigma0_conventional(idx_int_ISR(start_sigma0:finish_sigma0)));
    %----------------- SWH ------------------------------------
    ISR_L2_SWH_conventional_smoothed(i_wfm) = nanmean(ISR_L2_SWH_conventional(idx_int_ISR(start_SWH:finish_SWH)));
    % ---------- COR: pearson correlation coefficient ---------
    ISR_L2_COR_conventional_smoothed(i_wfm) = nanmean(ISR_L2_COR_conventional(idx_int_ISR(start_COR:finish_COR)));
    %-------------------- ACDC --------------------------------------------
    % compute the errors w.r.t fitting on the data using a smooth
    % function
    %------------------- SSH ----------------------------------
    ISR_L2_SSH_ACDC_smoothed(i_wfm)   = nanmean(ISR_L2_SSH_ACDC(idx_int_ISR(start_SSH:finish_SSH)));
    %----------------- sigma0 ---------------------------------
    ISR_L2_sigma0_ACDC_smoothed(i_wfm) = nanmean(ISR_L2_sigma0_ACDC(idx_int_ISR(start_sigma0:finish_sigma0)));
    %------------------- SWH ----------------------------------
    ISR_L2_SWH_ACDC_smoothed(i_wfm)   = nanmean(ISR_L2_SWH_ACDC(idx_int_ISR(start_SWH:finish_SWH)));
    % ---------- COR: pearson correlation coefficient ---------
    ISR_L2_COR_ACDC_smoothed(i_wfm) = nanmean(ISR_L2_COR_ACDC(idx_int_ISR(start_COR:finish_COR)));
    
end
   




%% --------------------------- SSH ----------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- --- --------------------------------------------$\n');
fprintf(fid,'$---------------- SSH --------------------------------------------$\n');
fprintf(fid,'$---------------- --- --------------------------------------------$\n');

if ~isempty(filename_L2_ESA)
    %remove outliers:
    ESA_L2_SSH_r1_filtered=ESA_L2_SSH_r1(ESA_indexes_int);
    %idx_outliers_corr=ESA_L2_SSH_r1_filtered>(nanmean(ESA_L2_SSH_r1_filtered)+threshold_std*nanstd(ESA_L2_SSH_r1_filtered)) | ESA_L2_SSH_r1_filtered<(nanmean(ESA_L2_SSH_r1_filtered)-threshold_std*nanstd(ESA_L2_SSH_r1_filtered));
    
    if ESA_num_surfaces_filtered==ISR_num_surfaces_filtered
        fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
        fprintf(fid,'$---------------- RETRACKERS COMPARISON ESA -------------------------------------------$\n');
        fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');        
        fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISR_L2_SSH_conventional_filtered=ISR_L2_SSH_conventional(ISR_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISR_L2_SSH_conventional;
        res.SSH.CONVENTIONAL.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISR_L2_SSH_conventional_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-ISR (CONVENTIONAL-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.CONVENTIONAL.RMSE_error_L2);
        res.SSH.CONVENTIONAL.mean_error_L2=(nanmean(abs(ESA_L2_SSH_r1_filtered(~idx_outliers)-ISR_L2_SSH_conventional_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-ISR (CONVENTIONAL-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.CONVENTIONAL.mean_error_L2);
        
        fprintf(fid,'$---------------- ACDC RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISR_L2_SSH_ACDC_filtered=ISR_L2_SSH_ACDC(ISR_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISR_L2_SSH_ACDC;
        res.SSH.ACDC.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISR_L2_SSH_ACDC_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-ISR (ACDC-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ACDC.RMSE_error_L2);
        res.SSH.ACDC.mean_error_L2=(nanmean(abs(ESA_L2_SSH_r1_filtered(~idx_outliers)-ISR_L2_SSH_ACDC_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-ISR (ACDC-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ACDC.mean_error_L2);
    end
end
fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS COMPARISON ACDC-CONVENTIONAL -----------------------------$\n');
fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
%----------------- CONVENTIONAL -------------------------------------------
ISR_L2_SSH_conventional_filtered=ISR_L2_SSH_conventional(ISR_indexes_int);
ISR_L2_SSH_ACDC_filtered=ISR_L2_SSH_ACDC(ISR_indexes_int);

idx_outliers=idx_outliers_ISR_L2_SSH_conventional | idx_outliers_ISR_L2_SSH_ACDC;
res.SSH.CONVENTIONAL.RMSE_error_ACDC=sqrt(nanmean((ISR_L2_SSH_conventional_filtered(~idx_outliers)-ISR_L2_SSH_ACDC_filtered(~idx_outliers)).^2));
fprintf(fid,'RMSE error SSH CONVENTIONAL-ACDC [m]: '); fprintf(fid,'%.18g\n',res.SSH.CONVENTIONAL.RMSE_error_ACDC);
res.SSH.CONVENTIONAL.mean_error_ACDC=(nanmean((ISR_L2_SSH_conventional_filtered(~idx_outliers)-ISR_L2_SSH_ACDC_filtered(~idx_outliers))));
fprintf(fid,'Mean error SSH CONVENTIONAL-ACDC [m]: '); fprintf(fid,'%.18g\n',res.SSH.CONVENTIONAL.mean_error_ACDC);
        

fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');
res.SSH.CONVENTIONAL.mean_error_fitting=nanmean(ISR_L2_SSH_conventional_filtered(~idx_outliers_ISR_L2_SSH_conventional)-ISR_L2_SSH_conventional_smoothed(~idx_outliers_ISR_L2_SSH_conventional));
fprintf(fid,'Mean error SSH fitting (CONVENTIONAL-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.CONVENTIONAL.mean_error_fitting);
res.SSH.CONVENTIONAL.rmse_fitting=sqrt(nanmean((ISR_L2_SSH_conventional_filtered(~idx_outliers_ISR_L2_SSH_conventional)-ISR_L2_SSH_conventional_smoothed(~idx_outliers_ISR_L2_SSH_conventional)).^2));
fprintf(fid,'RMSE SSH fitting (CONVENTIONAL-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.CONVENTIONAL.rmse_fitting);

fprintf(fid,'$---------------- ACDC RETRACKER -------------------------------------------$\n');
res.SSH.ACDC.mean_error_fitting=nanmean(ISR_L2_SSH_ACDC_filtered(~idx_outliers_ISR_L2_SSH_ACDC)-ISR_L2_SSH_ACDC_smoothed(~idx_outliers_ISR_L2_SSH_ACDC));
fprintf(fid,'Mean error SSH fitting (ACDC-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ACDC.mean_error_fitting);
res.SSH.ACDC.rmse_fitting=sqrt(nanmean((ISR_L2_SSH_ACDC_filtered(~idx_outliers_ISR_L2_SSH_ACDC)-ISR_L2_SSH_ACDC_smoothed(~idx_outliers_ISR_L2_SSH_ACDC)).^2));
fprintf(fid,'RMSE SSH fitting (ACDC-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ACDC.rmse_fitting);

if ~isempty(filename_L2_ESA)
    fprintf(fid,'$---------------- ESA RETRACKER -------------------------------------------$\n');
    res.SSH.ESA_L2.mean_error_fitting=nanmean(ESA_L2_SSH_r1_filtered(~idx_outliers_ESA_L2_SSH)-ESA_L2_SSH_r1_smoothed(~idx_outliers_ESA_L2_SSH));
    fprintf(fid,'Mean error SSH fitting (L2-ESA) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ESA_L2.mean_error_fitting);
    res.SSH.ESA_L2.rmse_fitting=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers_ESA_L2_SSH)-ESA_L2_SSH_r1_smoothed(~idx_outliers_ESA_L2_SSH)).^2));
    fprintf(fid,'RMSE SSH fitting (L2-ESA) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ESA_L2.rmse_fitting);
end


%---------------- Plot comparison ESA & ISR -------------------------------
figure;  
text_in_textbox={''};
legend_text={''};

%------------------- Conventional------------------------------------------
plot(ISR_lat_surf(ISR_indexes_int),ISR_L2_SSH_conventional(ISR_indexes_int),linestyle_conventional);
legend_text=[legend_text,strcat('ISR',{' '},title_name_conventional)];
text_in_textbox=[text_in_textbox,...
    strcat({'ISR '},title_name_conventional,' --> RMSE [m]:',{' '},num2str(res.SSH.CONVENTIONAL.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.CONVENTIONAL.mean_error_fitting,'%.4g'))];
hold on;
grid on;
%------------------ ACDC --------------------------------------------------
plot(ISR_lat_surf(ISR_indexes_int),ISR_L2_SSH_ACDC(ISR_indexes_int),linestyle_ACDC);
legend_text=[legend_text,strcat('ISR ',{' '},title_name_ACDC)];
text_in_textbox=[text_in_textbox,...
    strcat({'ISR '},title_name_ACDC,' --> RMSE [m]:',{' '},num2str(res.SSH.ACDC.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.ACDC.mean_error_fitting,'%.4g'))];

if ~isempty(filename_L2_ESA)
    plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1(ESA_indexes_int),linestyle_ESA);
    hold on;
    legend_text=[legend_text,{'ESA'}];
    text_in_textbox=[text_in_textbox,...
        strcat('ESA --> RMSE [m]:',{' '},num2str(res.SSH.ESA_L2.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.ESA_L2.mean_error_fitting,'%.4g'))];
    title(strcat('Comparison',{' '},upper(sh_name_nc),': ESA & isardSAT'));
else
    title(strcat('Comparison',{' '},upper(sh_name_nc),': isardSAT'));
end
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('Latitude [deg.]'); ylabel(strcat(upper(sh_name_nc),' [m]'));
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1)+0.015,0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
%print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISR),'.nc','_'),'Comparison_SSH_ESA_ISR.png'))
print('-dpng',strcat(output_path,name_file_L1B_ISR,'_cmp_SSH.png'))
    
%% ------------------------------ SIGMA0 ----------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$---------------- SIGMA0  --------------------------------------------$\n');
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
if ~isempty(filename_L2_ESA)
    ESA_L2_sigma0_r1_filtered=ESA_L2_sigma0_r1(ESA_indexes_int);
    %idx_outliers_ESA_L2_SSH=zeros(1,ESA_num_surfaces_filtered);
    if ESA_num_surfaces_filtered==ISR_num_surfaces_filtered
        fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
        fprintf(fid,'$---------------- RETRACKERS COMPARISON ESA -------------------------------------------$\n');
        fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
        %----------------- CONVENTIONAL -------------------------------------
        
        fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISR_L2_sigma0_conventional_filtered=ISR_L2_sigma0_conventional(ISR_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_sigma0 | idx_outliers_ISR_L2_sigma0_conventional;
        res.SIGMA0.CONVENTIONAL.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISR_L2_sigma0_conventional_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISR (CONVENTIONAL-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.CONVENTIONAL.RMSE_error_L2);
        res.SIGMA0.CONVENTIONAL.mean_error_L2=(nanmean(abs(ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISR_L2_sigma0_conventional_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISR (CONVENTIONAL-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.CONVENTIONAL.mean_error_L2);
        %----------------- ACDC -------------------------------------        
        fprintf(fid,'$---------------- ACDC RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISR_L2_sigma0_ACDC_filtered=ISR_L2_sigma0_ACDC(ISR_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_sigma0 | idx_outliers_ISR_L2_sigma0_ACDC;
        res.SIGMA0.ACDC.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISR_L2_sigma0_ACDC_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISR (ACDC-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ACDC.RMSE_error_L2);
        res.SIGMA0.ACDC.mean_error_L2=(nanmean(abs(ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISR_L2_sigma0_ACDC_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISR (ACDC-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ACDC.mean_error_L2);
    end
end
fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS COMPARISON ACDC-CONVENTIONAL -----------------------------$\n');
fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
%----------------- CONVENTIONAL -------------------------------------
fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');
%----------- Comparison with ESA L2 ---------------------------
ISR_L2_sigma0_conventional_filtered=ISR_L2_sigma0_conventional(ISR_indexes_int);
ISR_L2_sigma0_ACDC_filtered=ISR_L2_sigma0_ACDC(ISR_indexes_int);
idx_outliers=idx_outliers_ISR_L2_sigma0_conventional | idx_outliers_ISR_L2_sigma0_ACDC;
res.SIGMA0.CONVENTIONAL.RMSE_error_ACDC=sqrt(nanmean((ISR_L2_sigma0_conventional_filtered(~idx_outliers)-ISR_L2_sigma0_ACDC_filtered(~idx_outliers)).^2));
fprintf(fid,'RMSE error sigma0 CONVENTIONAL-ACDC [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.CONVENTIONAL.RMSE_error_ACDC);
res.SIGMA0.CONVENTIONAL.mean_error_ACDC=(nanmean(abs(ISR_L2_sigma0_conventional_filtered(~idx_outliers)-ISR_L2_sigma0_ACDC_filtered(~idx_outliers))));
fprintf(fid,'Mean error sigma0 CONVENTIONAL-ACDC [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.CONVENTIONAL.mean_error_ACDC);

fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');
res.SIGMA0.CONVENTIONAL.mean_error_fitting=nanmean(ISR_L2_sigma0_conventional(ISR_indexes_int)-ISR_L2_sigma0_conventional_smoothed);
fprintf(fid,'Mean error sigma0 fitting (CONVENTIONAL-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.CONVENTIONAL.mean_error_fitting);
res.SIGMA0.CONVENTIONAL.rmse_fitting=sqrt(nanmean((ISR_L2_sigma0_conventional(ISR_indexes_int)-ISR_L2_sigma0_conventional_smoothed).^2));
fprintf(fid,'RMSE sigma0 fitting (CONVENTIONAL-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.CONVENTIONAL.rmse_fitting);

fprintf(fid,'$---------------- ACDC RETRACKER -------------------------------------------$\n');
res.SIGMA0.ACDC.mean_error_fitting=nanmean(ISR_L2_sigma0_ACDC(ISR_indexes_int)-ISR_L2_sigma0_ACDC_smoothed);
fprintf(fid,'Mean error sigma0 fitting (ACDC-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ACDC.mean_error_fitting);
res.SIGMA0.ACDC.rmse_fitting=sqrt(nanmean((ISR_L2_sigma0_ACDC(ISR_indexes_int)-ISR_L2_sigma0_ACDC_smoothed).^2));
fprintf(fid,'RMSE sigma0 fitting (ACDC-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ACDC.rmse_fitting);

if ~isempty(filename_L2_ESA)
    fprintf(fid,'$---------------- ESA RETRACKER -------------------------------------------$\n');
    res.SIGMA0.ESA_L2.mean_error_fitting=nanmean(ESA_L2_sigma0_r1_filtered(~idx_outliers_ESA_L2_sigma0)-ESA_L2_sigma0_r1_smoothed(~idx_outliers_ESA_L2_sigma0));
    fprintf(fid,'Mean error sigma0 fitting (L2-ESA) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ESA_L2.mean_error_fitting);
    res.SIGMA0.ESA_L2.rmse_fitting=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers_ESA_L2_sigma0)-ESA_L2_sigma0_r1_smoothed(~idx_outliers_ESA_L2_sigma0)).^2));
    fprintf(fid,'RMSE sigma0 fitting (L2-ESA) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ESA_L2.rmse_fitting);
end


%---------------- Plot comparison ESA & ISR -------------------------------
figure;
text_in_textbox={''};
legend_text={''};

plot(ISR_lat_surf(ISR_indexes_int),ISR_L2_sigma0_conventional(ISR_indexes_int),linestyle_conventional);
legend_text=[legend_text,strcat('ISR',{' '},title_name_conventional)];
text_in_textbox=[text_in_textbox,...
    strcat({'ISR '},title_name_conventional,' --> RMSE [dB]:',{' '},num2str(res.SIGMA0.CONVENTIONAL.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.CONVENTIONAL.mean_error_fitting,'%.4g'))];
hold on;
grid on;

plot(ISR_lat_surf(ISR_indexes_int),ISR_L2_sigma0_ACDC(ISR_indexes_int),linestyle_ACDC);
legend_text=[legend_text,strcat('ISR',{' '},title_name_ACDC)];
text_in_textbox=[text_in_textbox,...
    strcat({'ISR '},title_name_ACDC,' --> RMSE [dB]:',{' '},num2str(res.SIGMA0.ACDC.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.ACDC.mean_error_fitting,'%.4g'))];

if ~isempty(filename_L2_ESA)
    plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),linestyle_ESA);
    hold on;
    legend_text=[legend_text,{'ESA'}];
    text_in_textbox=[text_in_textbox,...
        strcat('ESA --> RMSE [dB]:',{' '},num2str(res.SIGMA0.ESA_L2.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.ESA_L2.mean_error_fitting,'%.4g'))];
    title(strcat('Comparison \sigma^0: ESA & isardSAT'),'Interpreter','Tex');
else
    title(strcat('Comparison \sigma^0: isardSAT'),'Interpreter','Tex');
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]','Interpreter','Tex');
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1)+0.015,0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
print('-dpng',strcat(output_path,name_file_L1B_ISR,'_cmp_SIG0.png'))


%% ------------------------------ SWH -------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$----------------   SWH   --------------------------------------------$\n');
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
ISR_L2_SWH_conventional_filtered=ISR_L2_SWH_conventional(ISR_indexes_int);
ISR_L2_SWH_ACDC_filtered=ISR_L2_SWH_ACDC(ISR_indexes_int);

fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS COMPARISON CONVENTIONAL-ACDC -----------------------------$\n');
fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');
idx_outliers=idx_outliers_ISR_L2_SWH_conventional | idx_outliers_ISR_L2_SWH_ACDC;
res.SWH.CONVENTIONAL.RMSE_error_ACDC=sqrt(nanmean((ISR_L2_SWH_conventional_filtered(~idx_outliers)-ISR_L2_SWH_ACDC_filtered(~idx_outliers)).^2));
fprintf(fid,'RMSE error SWH CONVENTIONAL-ACDC [m]: '); fprintf(fid,'%.18g\n',res.SWH.CONVENTIONAL.RMSE_error_ACDC);
res.SWH.CONVENTIONAL.mean_error_ACDC=(nanmean(abs(ISR_L2_SWH_conventional_filtered(~idx_outliers)-ISR_L2_SWH_ACDC_filtered(~idx_outliers))));
fprintf(fid,'Mean error SWH CONVENTIONAL-ACDC [m]: '); fprintf(fid,'%.18g\n',res.SWH.CONVENTIONAL.mean_error_ACDC);

fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');

fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');
res.SWH.CONVENTIONAL.mean_error_fitting=nanmean(ISR_L2_SWH_conventional_filtered(~idx_outliers_ISR_L2_SWH_conventional)-ISR_L2_SWH_conventional_smoothed(~idx_outliers_ISR_L2_SWH_conventional));
fprintf(fid,'Mean error SWH fitting (CONVENTIONAL-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.CONVENTIONAL.mean_error_fitting);
res.SWH.CONVENTIONAL.rmse_fitting=sqrt(nanmean((ISR_L2_SWH_conventional_filtered(~idx_outliers_ISR_L2_SWH_conventional)-ISR_L2_SWH_conventional_smoothed(~idx_outliers_ISR_L2_SWH_conventional)).^2));
fprintf(fid,'RMSE SWH fitting (CONVENTIONAL-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.CONVENTIONAL.rmse_fitting);

fprintf(fid,'$------------------------ ACDC RETRACKER -------------------------------------------$\n');
res.SWH.ACDC.mean_error_fitting=nanmean(ISR_L2_SWH_ACDC_filtered(~idx_outliers_ISR_L2_SWH_ACDC)-ISR_L2_SWH_ACDC_smoothed(~idx_outliers_ISR_L2_SWH_ACDC));
fprintf(fid,'Mean error SWH fitting (ACDC-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ACDC.mean_error_fitting);
res.SWH.ACDC.rmse_fitting=sqrt(nanmean((ISR_L2_SWH_ACDC_filtered(~idx_outliers_ISR_L2_SWH_ACDC)-ISR_L2_SWH_ACDC_smoothed(~idx_outliers_ISR_L2_SWH_ACDC)).^2));
fprintf(fid,'RMSE SWH fitting (ACDC-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ACDC.rmse_fitting);


%----------------------------- ploting ------------------------------------
figure; 
legend_text={''};
text_in_textbox={''};

plot(ISR_lat_surf(ISR_indexes_int),abs(ISR_L2_SWH_conventional(ISR_indexes_int)),linestyle_conventional);
hold on;
grid on;
%plot(ISR_lat_surf(ISR_indexes_int),abs(ISR_L2_SWH_conventional_smoothed),linestyle_conventional_smooth);
legend_text=[legend_text,strcat('ISR',{' '},title_name_conventional)];%,strcat('L2-ISR',{' '},title_name_conventional,' (smoothed)')];
text_in_textbox=[text_in_textbox,...
    strcat({'ISR '},title_name_conventional,' --> RMSE [m]:',{' '},num2str(res.SWH.CONVENTIONAL.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SWH.CONVENTIONAL.mean_error_fitting,'%.4g'))];

plot(ISR_lat_surf(ISR_indexes_int),abs(ISR_L2_SWH_ACDC(ISR_indexes_int)),linestyle_ACDC);
%plot(ISR_lat_surf(ISR_indexes_int),abs(ISR_L2_SWH_ACDC_smoothed),linestyle_ACDC_smooth);
legend_text=[legend_text,strcat('ISR',{' '},title_name_ACDC)];%,strcat('ISR',{' '},title_name_conventional,' (smoothed)')];
text_in_textbox=[text_in_textbox,...
    strcat({'ISR '},title_name_ACDC,' --> RMSE [m]:',{' '},num2str(res.SWH.ACDC.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SWH.ACDC.mean_error_fitting,'%.4g'))];

title(strcat('Comparison SWH: isardSAT'),'Interpreter','Tex');
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('Latitude [deg.]'); ylabel('SWH [m]','Interpreter','Tex');
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1)+0.015,0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
%print('-dpng',strcat(output_path,strrep(char(name_file_L1B_ISR),'.nc','_'),'Comparison_ISR_SWH_smoothed.png'))
print('-dpng',strcat(output_path,name_file_L1B_ISR,'_ISR_SWH.png'))

%% ---------- COR:PEARSON CORRELATION COEFFICIENT -------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$--------   COR:PEARSON CORRELATION COEFFICIENT   --------------------$\n');
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
ISR_L2_COR_conventional_filtered=ISR_L2_COR_conventional(ISR_indexes_int);
ISR_L2_COR_ACDC_filtered=ISR_L2_COR_ACDC(ISR_indexes_int);

fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS COMPARISON CONVENTIONAL-ACDC -----------------------------$\n');
fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');

idx_outliers=idx_outliers_ISR_L2_COR_conventional | idx_outliers_ISR_L2_COR_ACDC;
res.COR.CONVENTIONAL.RMSE_error_ACDC=sqrt(nanmean((ISR_L2_COR_conventional_filtered(~idx_outliers)-ISR_L2_COR_ACDC_filtered(~idx_outliers)).^2));
fprintf(fid,'RMSE error COR CONVENTIONAL-ACDC [percentage]: '); fprintf(fid,'%.18g\n',res.COR.CONVENTIONAL.RMSE_error_ACDC);
res.COR.CONVENTIONAL.mean_error_ACDC=(nanmean(abs(ISR_L2_COR_conventional_filtered(~idx_outliers)-ISR_L2_COR_ACDC_filtered(~idx_outliers))));
fprintf(fid,'Mean error COR CONVENTIONAL-ACDC [percentage]: '); fprintf(fid,'%.18g\n',res.COR.CONVENTIONAL.mean_error_ACDC);

fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
fprintf(fid,'$---------------- CONVENTIONAL RETRACKER -------------------------------------------$\n');
res.COR.CONVENTIONAL.mean=nanmean(ISR_L2_COR_conventional(ISR_indexes_int));
fprintf(fid,'Mean value Pearson Correlation Coefficient (CONVENTIONAL-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.CONVENTIONAL.mean);
res.COR.CONVENTIONAL.mean_error_fitting=nanmean(ISR_L2_COR_conventional(ISR_indexes_int)-ISR_L2_COR_conventional_smoothed);
fprintf(fid,'Mean error Pearson Correlation Coefficient fitting (CONVENTIONAL-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.CONVENTIONAL.mean_error_fitting);
res.COR.CONVENTIONAL.rmse_fitting=sqrt(nanmean((ISR_L2_COR_conventional(ISR_indexes_int)-ISR_L2_COR_conventional_smoothed).^2));
fprintf(fid,'RMSE Pearson Correlation Coefficient fitting (CONVENTIONAL-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.CONVENTIONAL.rmse_fitting);
fprintf(fid,'$---------------- ACDC RETRACKER -------------------------------------------$\n');
res.COR.ACDC.mean=nanmean(ISR_L2_COR_ACDC(ISR_indexes_int));
fprintf(fid,'Mean value Pearson Correlation Coefficient (ACDC-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ACDC.mean);
res.COR.ACDC.mean_error_fitting=nanmean(ISR_L2_COR_ACDC(ISR_indexes_int)-ISR_L2_COR_ACDC_smoothed);
fprintf(fid,'Mean error Pearson Correlation Coefficient fitting (ACDC-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ACDC.mean_error_fitting);
res.COR.ACDC.rmse_fitting=sqrt(nanmean((ISR_L2_COR_ACDC(ISR_indexes_int)-ISR_L2_COR_ACDC_smoothed).^2));
fprintf(fid,'RMSE Pearson Correlation Coefficient fitting (ACDC-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ACDC.rmse_fitting);
%----------------------------- ploting ------------------------------------
figure; 
legend_text={''};
text_in_textbox={''};

plot(ISR_lat_surf(ISR_indexes_int),abs(ISR_L2_COR_conventional(ISR_indexes_int)),linestyle_conventional);
hold on;
grid on;
%plot(ISR_lat_surf(ISR_indexes_int),abs(ISR_L2_COR_conventional_smoothed),linestyle_conventional_smooth);
legend_text=[legend_text,strcat('ISR',{' \rho_{pearson} '},title_name_conventional)];%,strcat('ISR',{' \rho_{pearson} '},title_name_conventional,' (smoothed)')];
text_in_textbox=[text_in_textbox,...
    strcat({'ISR '},title_name_conventional,' --> RMSE [%]:',{' '},num2str(res.COR.CONVENTIONAL.rmse_fitting,'%.4g'))];%,', Mean [%]:',{' '},num2str(res.COR.CONVENTIONAL.mean_error_fitting,'%.4g'))];

plot(ISR_lat_surf(ISR_indexes_int),abs(ISR_L2_COR_ACDC(ISR_indexes_int)),linestyle_ACDC);
%plot(ISR_lat_surf(ISR_indexes_int),abs(ISR_L2_COR_ACDC_smoothed),linestyle_ACDC_smooth);
legend_text=[legend_text,strcat('ISR',{' \rho_{pearson} '},title_name_ACDC)];%,strcat('ISR',{' \rho_{pearson} '},title_name_ACDC,' (smoothed)')];
text_in_textbox=[text_in_textbox,...
    strcat({'ISR '},title_name_ACDC,' --> RMSE [%]:',{' '},num2str(res.COR.ACDC.rmse_fitting,'%.4g'))];%,', Mean [%]:',{' '},num2str(res.COR.ACDC.mean_error_fitting,'%.4g'))];

title(strcat('Comparison goodness of fitting (gof): isardSAT'),'Interpreter','Tex');
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('Latitude [deg.]'); ylabel('gof [%]','Interpreter','Tex');
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1)+0.015,0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
%print('-dpng',strcat(output_path,strrep(char(name_file_L1B_ISR),'.nc','_'),'Comparison_ISR_COR_RHO_smoothed.png'))
print('-dpng',strcat(output_path,name_file_L1B_ISR,'_ISR_COR.png'))

%% --------------------------- Saving Information -------------------------
fclose(fid);
close all;
%save(strcat(output_path,strrep(char(name_file_L1B_ISR),'.nc','_'),'L2_Evaluation.mat'),'res');
save(strcat(output_path,name_file_L1B_ISR,'_L2_Evaluation.mat'),'res');

end