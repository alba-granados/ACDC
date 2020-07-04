% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT Ltd. 
% --------------------------------------------------------
%
% Ground Processor Prototype for altimetry missions:
% CryoSat-2 SARIn
% CryoSat-2 SAR
% Sentinel 3 SAR
% Sentinel 6 RAW
% Sentinel 6 RMC
% 
% ---------------------------------------------------------
% Inputs: 
% L0 (ISP)+ orbit file + attitude file + meteo file
% L1A (FBR in case of CryoSat-2)
% L1B-S
% L1B
% 
% Output
%   NetCDF L1A, L1B-S, L1B and/or L2
%
% ----------------------------------------------------------
% 
% Authors:  Eduard makhoul / isardSAT
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling GPP_bulk_processing_paralelization('CR2','./inputs/', './results/',1,0)
% Missions accepted: 'CR2' , 'S3'. 'S6' to be added
% 1.0
% 1.1 Changed FBR for L1A in the name of vars. Added S3 elseif for reading files
% 1.2 Added error log whe L1B processing crashes
% 1.3 options.writing_flags = [writing_L1B_pLRM writing_L1BS writing_L1B]
% and options.plotting_flags = [Stacks_flag Waveforms_flag] 
% 1.4 fid_log added to the files structure to write info in L1B chain
function GPP_bulk_processing_paralelization(mission, input_path,output_path,num_pools,targz_option_active, options)
clear global
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');
if(nargin < 6)
    options.writting_flag = [0 0 0 0]; % L1BS, L1B, pLRM, KML
    options.plotting_flag = [0 0 0]; % Stacks, L1B waveforms, track
    options.axes = [];
    options.wd_axes = [];
    options.GUI_flag = 0;
end
     
version_matlab=version;
ttotal=tic;
log_flag = 1;

if(log_flag)
    filesBulk.fid_log = fopen([output_path 'LogError.txt'],'a+');
end

%% ----- Include the different subfolders in the search path --------------
FolderInfo=dir;
FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
aux=struct2cell(FolderInfo); %into a cell array where first row is 
folders=(aux(1,[FolderInfo.isdir]))'; %name is the first row
clear aux;
folders=strcat('./',folders(~strcmp(folders,['.'])&~strcmp(folders,['..'])&~strcmp(folders,['.svn'])&~strcmp(folders,['inputs'])));
for i_folder=1:length(folders)
    addpath(genpath(char(folders(i_folder))));
end
if(sum(options.plotting_flag)>0)
    set_default_plot;
end
%% --------- Define Paths -------------------------------------------------
%fixing '/' '\'
if(strfind(input_path, '/'))
    if(~strcmp(input_path(end), '/'))
        input_path = [input_path '/'];
    end
elseif(strfind( input_path, '\'))
    if(~strcmp(input_path(end), '\'))
        input_path = [input_path '\'];
    end
end
if(strfind(output_path, '/'))
    if(~strcmp(output_path(end), '/'))
        output_path = [output_path '/'];
    end
elseif(strfind( output_path, '\'))
    if(~strcmp(output_path(end), '\'))
        output_path = [output_path '\'];
    end
end
filesBulk.inputPath       =   input_path;
disp(strcat('Input_path: ',filesBulk.inputPath))
filesBulk.resultPath      =   output_path;
disp(strcat('Result_path: ',filesBulk.resultPath))
mkdir(filesBulk.resultPath);
mkdir([filesBulk.resultPath 'data/']);
mkdir([filesBulk.resultPath 'plots/']);
filesBulk.inputFiles      =   dir(filesBulk.inputPath);
filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles          =   length(filesBulk.indexFiles);             % number of input files

aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
if(strcmp(mission,'CR2'))
    if targz_option_active 
       filesBulk.filterDATAFILES=(~cellfun(@isempty,strfind(aux,'TGZ')));
       filesBulk.indexFilesL1A=find(filesBulk.filterDATAFILES);
    else
       filesBulk.indexFilesL1A=find(~cellfun(@isempty,strfind(aux,'DBL')));
       filesBulk.filterDATAFILES=(~cellfun(@isempty,strfind(aux,'DBL'))) | (~cellfun(@isempty,strfind(aux,'HDR')));
    end
elseif(strcmp(mission,'S3'))
    if targz_option_active 
       filesBulk.filterDATAFILES=(~cellfun(@isempty,strfind(aux,'TGZ')));
       filesBulk.indexFilesL1A=find(filesBulk.filterDATAFILES);
    else
       filesBulk.indexFilesL1A=find(~cellfun(@isempty,strfind(aux,'_l1a.nc')));
       filesBulk.filterDATAFILES=(~cellfun(@isempty,strfind(aux,'_l1a.nc'))) | (~cellfun(@isempty,strfind(aux,'HDR')));
    end
end
filesBulk.nFilesL1A=length(filesBulk.indexFilesL1A);
filesBulk.L1AFiles=filesBulk.inputFiles(filesBulk.indexFilesL1A);

disp('Total number of L1A to be processed');
disp(num2str(filesBulk.nFilesL1A))

% %---- Display the main processing configuration parameters ------
% CNF_file=[filesBulk.inputPath
% filesBulk.inputFiles(~cellfun(@isempty,strfind(aux,'cnf_file'))).name]; %
% alba: these paths are not correct. See ../config/ in
% L1B_processing_ACDC_bis.m
% CHD_file=[filesBulk.inputPath filesBulk.inputFiles(~cellfun(@isempty,strfind(aux,'chd_file'))).name];
% run(CHD_file);
% run(CNF_file);
% global hamming_window_cnf height_rate_application_cnf FAI_application_cnf CAL2_flag_cnf CAL1p2p_flag_cnf zp_fact_range_cnf
% global apply_stack_mask_cnf avoid_beams_mask_allzeros_cnf avoid_noisy_beams_cnf method_noise_est_cnf param_std_noise_thres_cnf noise_estimation_window
% disp('------------ Main processing configuration parameters -------------')
% disp(strcat('height_rate_application_cnf: ',num2str(height_rate_application_cnf)));
% disp(strcat('FAI_application_cnf: ',num2str(FAI_application_cnf)));
% disp(strcat('CAL2_flag_cnf: ',num2str(CAL2_flag_cnf)));
% disp(strcat('CAL1p2p_flag_cnf: ',num2str(CAL1p2p_flag_cnf)));
% disp(strcat('hamming_window_cnf: ',num2str(hamming_window_cnf)));
% disp(strcat('zp_fact_range_cnf: ',num2str(zp_fact_range_cnf)));
% disp(strcat('apply_stack_mask_cnf: ',num2str(apply_stack_mask_cnf)));
% disp(strcat('avoid_beams_mask_allzeros_cnf: ',num2str(avoid_beams_mask_allzeros_cnf)));
% disp(strcat('avoid_noisy_beams_cnf: ',num2str(avoid_noisy_beams_cnf)));
% if avoid_noisy_beams_cnf    
%     disp(strcat('method_noise_est_cnf: ',method_noise_est_cnf)); 
%     disp(strcat('param_std_noise_thres_cnf: ',num2str(param_std_noise_thres_cnf)));   
%     disp(strcat('noise_estimation_window: ',num2str(noise_estimation_window)));   
% end
% disp('-------------------------------------------------------------------')

%% --------------- Run parallel processing --------------------------------
if num_pools~=1
    %create pools
    if str2double(version_matlab(end-5:end-2))>2013
        parpool(num_pools);    
    else
        matlabpool('open',num_pools);
    end
    %% ------------- Loop per file to be processed ------------------------
    parfor i_fileL1A_input=1:filesBulk.nFilesL1A
        try
            L1B_processing_ACDC_bis(filesBulk,i_fileL1A_input,targz_option_active,options);
            Logtime = datestr(now, 'yyyymmddTHHMMSS');
            fprintf(filesBulk.fid_log,'%s -> OKKKK!: File %s succesfully processed\n',Logtime, filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name);
        catch err
            if(log_flag)
                Logtime = datestr(now, 'yyyymmddTHHMMSS');
                fprintf(filesBulk.fid_log,'%s -> ERROR!: processing file %s\n',Logtime, filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name);
                fprintf(filesBulk.fid_log, '%s\n', err.getReport('extended', 'hyperlinks','off'));
            end
            continue;
        end
    end
    %close pools
    if str2double(version_matlab(end-5:end-2))>2013
        poolobj = gcp('nocreate');
        delete(poolobj);
    else
        matlabpool('close');
    end    
else
    for i_fileL1A_input=1:filesBulk.nFilesL1A
        try
            L1B_processing_ACDC_bis(filesBulk,i_fileL1A_input,targz_option_active,options);
            Logtime = datestr(now, 'yyyymmddTHHMMSS');
            fprintf(filesBulk.fid_log,'%s -> OKKKK!: File %s succesfully processed\n',Logtime, filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name);
         catch err
            if(log_flag)
                Logtime = datestr(now, 'yyyymmddTHHMMSS');
                fprintf(filesBulk.fid_log,'%s -> ERROR!: processing file %s\n',Logtime, filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name);
                fprintf(filesBulk.fid_log, '%s\n', err.getReport('extended', 'hyperlinks','off'));

            end
            continue;
        end
    end
end

time = toc(ttotal);
hours_processing = floor(time/3600);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time: ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);
if(log_flag)

    fprintf(filesBulk.fid_log,'Total processing time: %s hours %s minutes and %s seconds',num2str(hours_processing), num2str(minutes_processing),num2str(secs_processing));
    fclose (filesBulk.fid_log);
end

end
    
    
    
    


