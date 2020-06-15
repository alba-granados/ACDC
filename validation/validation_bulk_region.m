function L2_val = validation_bulk_region(input_path_L1B_ISD,input_path_L1A_ESA,input_path_L1B_ESA,input_path_L2_ESA,output_path_L2_ISD,path_comparison_results,mission,mode,RETRACKER)

surface_aligned=0;
retrack_flag=1;

filesBulk.inputPath       =   input_path_L1B_ISD;

% mkdir(path_comparison_results);

filesBulk.inputFiles      =   dir(filesBulk.inputPath);
filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles          =   length(filesBulk.indexFiles);             % number of input files
aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
filesBulk.indexFilesNC=find(~cellfun(@isempty,strfind(aux,'.nc')));
filesBulk.nFilesNC=length(filesBulk.indexFilesNC);
filesBulk.NCFiles=filesBulk.inputFiles(filesBulk.indexFilesNC);

i_files_valid = 1;
for i_file=1:filesBulk.nFilesNC
    filename_L1B_ISD    = char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name);
    data_string         = filename_L1B_ISD(17:17+30);
    filename_L1B_ISD    = strcat(input_path_L1B_ISD,filename_L1B_ISD);
    
    inputL1AESAFile    = dir(fullfile(input_path_L1A_ESA,['*' data_string(1:15) '*.DBL']));
    inputL1BESAFile    = dir(fullfile(input_path_L1B_ESA,['*' data_string(1:15) '*.DBL']));
    inputL2ESAFile     = dir(fullfile(input_path_L2_ESA,['*' data_string(1:13) '*.DBL']));
    
    if ~isempty(inputL1BESAFile) && ~isempty(inputL2ESAFile)
        indexL1AESA         = not([inputL1AESAFile.isdir]);
        filename_L1A_ESA    = strcat(input_path_L1A_ESA,char(inputL1AESAFile(indexL1AESA).name));
        SPH                 = readSPH(filename_L1A_ESA);
        
        indexL1BESA         = not([inputL1BESAFile.isdir]);
        filename_L1B_ESA    = strcat(input_path_L1B_ESA,char(inputL1BESAFile(indexL1BESA).name));
        disp(char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name))
        
        indexL2ESA          = not([inputL2ESAFile.isdir]);
        filename_L2_ESA     = strcat(input_path_L2_ESA,char(inputL2ESAFile(indexL2ESA).name));
        %filename_L1B_ESA=strcat('CS_LTA__SIR_SAR_1B_',data_string,'_C001.DBL');
        %filename_L1B_ESA=strcat(input_path_L1B_ESA,filename_L1B_ESA);
        %filename_L2_ESA=strcat('CS_LTA__SIR_SAR_2__',data_string,'_C001.DBL');
        %filename_L2_ESA=strcat(input_path_L2_ESA,filename_L2_ESA);
        
        [validation_structure(i_files_valid),L2(i_files_valid)] = L1B_validation (filename_L1B_ISD,filename_L1B_ESA,filename_L2_ESA,path_comparison_results,surface_aligned,retrack_flag,RETRACKER);
        
        L2_met(i_files_valid) = L2_metrics (L2(i_files_valid));
                
        %[L2_val(i_files_valid)] = create_NetCDF_L2_S3 (L2_met(i_files_valid),SPH,mission,mode,output_path_L2_ISD);
        
        %prepare_NetCDF_L2_S3 (L2_val(i_files_valid));

        %load(strcat(path_comparison_results,strrep(char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name),'.nc','_Evaluation.mat')));
        %validation_structure(i_files_valid)=res;
        %clear res;
        
% % %         GEO(i_files_valid)=validation_structure(i_files_valid).GEO;
% % %         ATT(i_files_valid)=validation_structure(i_files_valid).ATT;
% % %         ALT(i_files_valid)=validation_structure(i_files_valid).ALT;
% % %         SURF(i_files_valid)=validation_structure(i_files_valid).SURF;
        WINDELAY(i_files_valid)         = validation_structure(i_files_valid).WINDELAY;
        RETRACKED_RANGE(i_files_valid)  = validation_structure(i_files_valid).RETRACKED_RANGE.peak_detector;
        SIGMA0(i_files_valid)           = validation_structure(i_files_valid).SIGMA0.peak_detector;
        SSH(i_files_valid)              = validation_structure(i_files_valid).SSH.peak_detector;
        WVFMS(i_files_valid)            = validation_structure(i_files_valid).WVFMS;
        i_files_valid                   = i_files_valid+1;
        save(strcat(path_comparison_results,'Bulk_validation_information.mat'),'WINDELAY','RETRACKED_RANGE','SIGMA0','SSH','WVFMS');
%         save(strcat(path_comparison_results,'Bulk_validation_information.mat'),'GEO','ATT','ALT','SURF','WINDELAY','RETRACKED_RANGE','SIGMA0','SSH','WVFMS');
    end
end


%% ----------  Ploting ----------------------------------------------------

% retracked range
figure;
plot([RETRACKED_RANGE(:).RMSE_error],'-r'); 
title('RMSE error on retracked range L1B waveforms: ESA & ISD')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'RMSE_error_retrackedrange_ESA_ISD.png'))

figure;
plot([RETRACKED_RANGE(:).mean_error],'-r'); 
title('Mean error on retracked range L1B waveforms: ESA & ISD')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'Mean_error_retrackedrange_ESA_ISD.png'))

% SSH
figure;
plot([SSH(:).RMSE_error_L1B],'-r');
hold on; plot([SSH(:).RMSE_error_L2_nocorr],'-b');
title('RMSE error on SSH')
legend('L1B-ISD vs L1B-ESA (threshold retracker)','L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'RMSE_error_SSH_ESA_ISD.png'))

figure;
plot([SSH(:).mean_error_L1B],'-r');
hold on; plot([SSH(:).mean_error_L2_nocorr],'-b');
title('Mean error on SSH')
legend('L1B-ISD vs L1B-ESA (threshold retracker)','L1B-ISD (threshold retracker) vs L2-ESA')
xlabel('Track'); ylabel('[m]');
print('-dpng',strcat(path_comparison_results,'Mean_error_SSH_ESA_ISD.png'))

% Sigma0
figure;
plot([SIGMA0(:).RMSE_error_L2],'-r');
title('RMSE error on sigma0')
xlabel('Track'); ylabel('[dB]');
print('-dpng',strcat(path_comparison_results,'RMSE_error_sigma0_ESA_ISD.png'))

figure;
plot([SIGMA0(:).RMSE_error_L2],'-r');
title('RMSE error on sigma0')
xlabel('Track'); ylabel('[dB]');
print('-dpng',strcat(path_comparison_results,'Mean_error_sigma0_ESA_ISD.png'))

% waveforms
figure;
plot([WVFMS(:).RMSE_error_peak],'-r'); 
title('RMSE error on peak power L1B waveforms: ESA & ISD')
xlabel('Track'); ylabel('[dBW]');
print('-dpng',strcat(path_comparison_results,'RMSE_error_peakpower_ESA_ISD.png'))

figure;
plot([WVFMS(:).mean_error_peak],'-r'); 
title('Mean error on peak power L1B waveforms: ESA & ISD')
xlabel('Track'); ylabel('[dBW]');
print('-dpng',strcat(path_comparison_results,'Mean_error_peakpower_ESA_ISD.png'))



