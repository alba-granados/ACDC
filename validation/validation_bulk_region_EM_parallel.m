function validation_bulk_region_EM_parallel(input_path_L1B_ISD,input_path_L1B_ESA,input_path_L2_ESA,path_comparison_results,varargin)

version_matlab=version;
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<4 || nargin>6)
    error('Wrong number of input parameters');   
end
p = inputParser;
p.addParamValue('num_pools',1);
p.addParamValue('input_path_L2_STL',{''},@(x)ischar(x));
p.parse(varargin{:});
input_path_L2_STL=char(p.Results.input_path_L2_STL);
num_pools=(p.Results.num_pools);
clear p;

figures_visible=0;
surface_aligned=0;
retrack_flag=1;

if figures_visible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

filesBulk.inputPath       =   input_path_L1B_ISD;

mkdir(path_comparison_results);

filesBulk.inputFiles      =   dir(filesBulk.inputPath);
filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles          =   length(filesBulk.indexFiles);             % number of input files
aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
filesBulk.indexFilesNC=find(~cellfun(@isempty,strfind(aux,'.nc')));
filesBulk.nFilesNC=length(filesBulk.indexFilesNC);
filesBulk.NCFiles=filesBulk.inputFiles(filesBulk.indexFilesNC);

i_files_valid=1;
%% -------------- check available files -----------------------------------
%--------------------------------------------------------------------------
for i_file=1:filesBulk.nFilesNC
    filename_L1B_ISD=char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name);
    data_string=filename_L1B_ISD(17:17+30);
    filename_L1B_ISD=strcat(input_path_L1B_ISD,filename_L1B_ISD);       
    
    inputL1BESAFiles   = dir(fullfile(input_path_L1B_ESA,['*' data_string(1:15) '*.DBL']));
    inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' data_string(1:15) '*.DBL']));
    if isempty(input_path_L2_STL)
        inputL2STLFiles   = 1;
    else
        inputL2STLFiles   = dir(fullfile(input_path_L2_STL,['*' data_string(1:15) '*.nc']));
    end
    
    if ~isempty(inputL1BESAFiles) && ~isempty(inputL2ESAFiles) && ~isempty(inputL2STLFiles)
        filesBulk.indexFilesNC_valid(i_files_valid)=filesBulk.indexFilesNC(i_file);
        i_files_valid=i_files_valid+1;
    end
end

%% ------------------ RUN the validation for each available file ----------
%--------------------------------------------------------------------------
filesBulk.nFilesNC_valid=length(filesBulk.indexFilesNC_valid);
disp(strcat('Total number of valid files for evaluation: ',num2str(filesBulk.nFilesNC_valid)));
if num_pools~=1
    %create pools
    if str2double(version_matlab(end-5:end-2))>2013
        parpool(num_pools);    
    else
        matlabpool('open',num_pools);
    end
    %% ------------- Loop per file to be processed ------------------------
    parfor i_files_valid=1:filesBulk.nFilesNC_valid        
        run_L1B_validation(filesBulk,i_files_valid,input_path_L1B_ISD,...
            input_path_L1B_ESA,input_path_L2_ESA,...
            path_comparison_results,surface_aligned,retrack_flag,...
            'input_path_L2_STL',input_path_L2_STL,'figures_visible',figures_visible)              
    end
    %close pools
    if str2double(version_matlab(end-5:end-2))>2013
        poolobj = gcp('nocreate');
        delete(poolobj);
    else
        matlabpool('close');
    end    
else
    for i_files_valid=1:filesBulk.nFilesNC_valid
        run_L1B_validation(filesBulk,i_files_valid,input_path_L1B_ISD,...
            input_path_L1B_ESA,input_path_L2_ESA,...
            path_comparison_results,surface_aligned,retrack_flag,...
            'input_path_L2_STL',input_path_L2_STL,'figures_visible',figures_visible);   
    end
end

%% ------------- LOAD the output data structure ---------------------------
%--------------------------------------------------------------------------
%loading and reordering the data into a single array of structures
filesBulk.inputFilesEvaluation      =   dir(fullfile(path_comparison_results,'*_Evaluation.mat'));
for i_files_valid=1:filesBulk.nFilesNC_valid
    load(strcat(path_comparison_results,char(filesBulk.inputFilesEvaluation(i_files_valid).name)))
    GEO(i_files_valid)=res.GEO;
    ATT(i_files_valid)=res.ATT;
    ALT(i_files_valid)=res.ALT;
    SURF(i_files_valid)=res.SURF;
    WINDELAY(i_files_valid)=res.WINDELAY;
    RETRACKED_RANGE(i_files_valid)=res.RETRACKED_RANGE.peak_detector;
    SIGMA0(i_files_valid)=res.SIGMA0;
    SSH(i_files_valid)=res.SSH;
    WVFMS(i_files_valid)=res.WVFMS;
end
save(strcat(path_comparison_results,'Bulk_validation_information.mat'),'GEO','ATT','ALT','SURF','WINDELAY','RETRACKED_RANGE','SIGMA0','SSH','WVFMS');


%% ----------  Ploting ----------------------------------------------------
%% ---------------------------  SSH ---------------------------------------
%--------------------------------------------------------------------------
%-------------- Inter-comparison retrackers -------------------------------
peak_detector=[SSH(:).peak_detector];
figure;
plot([peak_detector(:).RMSE_error_L1B],'*r');
hold on; plot([peak_detector(:).RMSE_error_L2_nocorr],'ob');
title('RMSE error on SSH')
legend('L1B-ISD vs L1B-ESA (threshold retracker)','L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'RMSE_SSH_ESA_ISD.png'))

figure;
plot([peak_detector(:).mean_error_L1B],'*r');
hold on; plot([peak_detector(:).mean_error_L2_nocorr],'ob');
title('Mean error on SSH')
legend('L1B-ISD vs L1B-ESA (threshold retracker)','L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'Mean_error_SSH_ESA_ISD.png'))

%-------------- Fitting ---------------------------------------------------
ISD_peak_retracker=[peak_detector(:).ISD];
ESA_peak_retracker=[peak_detector(:).ESA];
ESA_retracker=[SSH(:).ESA_L2];
figure;
plot([ISD_peak_retracker(:).rmse_fitting],'ob');
hold on; plot([ESA_peak_retracker(:).rmse_fitting],'-m');
plot([ESA_retracker(:).rmse_fitting],'^r');
title('RMSE error on fitted SSH')
legend('L1B-ISD (threshold retracker)','L1B-ESA (threshold retracker)', 'L2-ESA')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'RMSE_fitted_SSH.png'))

figure;
plot([ISD_peak_retracker(:).mean_error_fitting],'ob');
hold on; plot([ESA_peak_retracker(:).mean_error_fitting],'-m');
plot([ESA_retracker(:).mean_error_fitting],'^r');
title('Mean error on fitted SSH')
legend('L1B-ISD (threshold retracker)','L1B-ESA (threshold retracker)', 'L2-ESA')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'Mean_error_fitted_SSH.png'))


%% ---------------------------  sigma0 ------------------------------------
%--------------------------------------------------------------------------
%-------------- Inter-comparison retrackers -------------------------------
peak_detector=[SIGMA0(:).peak_detector];
figure;
plot([peak_detector(:).RMSE_error_L2],'*r');
title('RMSE error on \sigma^0','Interpreter','Tex')
legend('L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[dB]');
print('-dpng',strcat(path_comparison_results,'RMSE_sigma0_ESA_ISD.png'))

figure;
plot([peak_detector(:).mean_error_L2],'*r');
title('Mean error on \sigma^0','Interpreter','Tex')
legend('L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[dB]');
print('-dpng',strcat(path_comparison_results,'Mean_error_sigma0_ESA_ISD.png'))

%-------------- Fitting ---------------------------------------------------
ISD_peak_retracker=[peak_detector(:).ISD];
ESA_retracker=[SIGMA0(:).ESA_L2];
figure;
plot([ISD_peak_retracker(:).rmse_fitting],'ob');
hold on;
plot([ESA_retracker(:).rmse_fitting],'^r');
title('RMSE error on fitted \sigma^0','Interpreter','Tex')
legend('L1B-ISD (threshold retracker)', 'L2-ESA')
xlabel('Track'); ylabel('[dB]');
print('-dpng',strcat(path_comparison_results,'RMSE_fitted_sigma0.png'))

figure;
plot([ISD_peak_retracker(:).mean_error_fitting],'ob');
hold on;
plot([ESA_retracker(:).mean_error_fitting],'^r');
title('Mean error on fitted \sigma^0','Interpreter','Tex')
legend('L1B-ISD (threshold retracker)', 'L2-ESA')
xlabel('Track'); ylabel('[dB]');
print('-dpng',strcat(path_comparison_results,'Mean_error_fitted_sigma0.png'))

%% ---------------------------  Peak power --------------------------------
%--------------------------------------------------------------------------
%-------------- Inter-comparison retrackers -------------------------------
figure;
plot([WVFMS(:).RMSE_error_peak],'*r');
title('RMSE error on peak power')
legend('L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[dB]');
print('-dpng',strcat(path_comparison_results,'RMSE_peakpower_ESA_ISD.png'))

figure;
plot([WVFMS(:).mean_error_peak],'*r');
title('Mean error on peak power')
legend('L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[dB]');
print('-dpng',strcat(path_comparison_results,'Mean_error_peakpower_ESA_ISD.png'))

%% ---------------------------  GEOLOCATION -------------------------------
%--------------------------------------------------------------------------
%-------------- Inter-comparison retrackers -------------------------------
GEO_LAT=[GEO(:).LAT];
GEO_LON=[GEO(:).LON];

figure;
plot([GEO(:).mean_error],'*r');
title('RMSE error on geolocation (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'Mean_geolocation_error_ESA_ISD.png'))

figure;
plot([GEO_LAT(:).RMSE_error],'*r');
title('RMSE error on latitude (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'RMSE_latitude_ESA_ISD.png'))

figure;
plot([GEO_LAT(:).mean_error],'*r');
title('Mean error on latitude (L1B-ESA vs L1B-ISD)')
legend('L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'Mean_error_latitude_ESA_ISD.png'))

figure;
plot([GEO_LON(:).RMSE_error],'*r');
title('RMSE error on longitude (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'RMSE_longitude_ESA_ISD.png'))

figure;
plot([GEO_LON(:).mean_error],'*r');
title('Mean error on longitude (L1B-ESA vs L1B-ISD)')
legend('L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'Mean_error_longitude_ESA_ISD.png'))

%% ---------------------  ATTITUDE INFORMATION ----------------------------
%--------------------------------------------------------------------------
%-------------- Inter-comparison retrackers -------------------------------
ATT_pitch=[ATT(:).pitch];
ATT_roll=[ATT(:).roll];
ATT_yaw=[ATT(:).yaw];
figure;
plot([ATT_pitch(:).RMSE_error],'*r');
title('RMSE error on pitch (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'RMSE_pitch_ESA_ISD.png'))

figure;
plot([ATT_pitch(:).mean_error],'*r');
title('Mean error on pitch (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'Mean_error_pitch_ESA_ISD.png'))

figure;
plot([ATT_roll(:).RMSE_error],'*r');
title('RMSE error on roll (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'RMSE_roll_ESA_ISD.png'))

figure;
plot([ATT_roll(:).mean_error],'*r');
title('Mean error on roll (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'Mean_error_roll_ESA_ISD.png'))

figure;
plot([ATT_yaw(:).RMSE_error],'*r');
title('RMSE error on yaw (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'RMSE_yaw_ESA_ISD.png'))

figure;
plot([ATT_yaw(:).mean_error],'*r');
title('Mean error on yaw (L1B-ESA vs L1B-ISD)')
xlabel('Track'); ylabel('[deg]');
print('-dpng',strcat(path_comparison_results,'Mean_error_yaw_ESA_ISD.png'))


