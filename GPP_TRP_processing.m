% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT Ltd.
% --------------------------------------------------------
%
% TRANSPONDER PROCESSOR S3 MPC:
%
% ---------------------------------------------------------
% Inputs:
% L1A
% ATMR file For tropo and iono corrections
% L2 data for the other geo corrections

%
% Output
%   .csv with TRP results
%
% ----------------------------------------------------------
%
% Authors:  Albert Garcia-Mondejar  / isardSAT

%%ENSURE MASK FLAG IN cnf_file.m IS SWITCHED BACK ON

function GPP_TRP_processing(mission, input_path, output_path, options)
clear global
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');
if(nargin < 4)
    options.writting_flag = [0 0 0 1]; % L1BS, L1B, pLRM, KML
    options.plotting_flag = [0 0 0]; % Stacks, L1B waveforms, track
    options.axes = [];
    options.wd_axes = [];
    options.GUI_flag = 0;
end



version_matlab=version;
ttotal=tic;
log_flag = 1;

global corr_mode mode_L1 output_File


corr_mode = 3; %[1,2,3,4] [pos,average,min,max]

pass_folder=dir('./201*');
indexaDirs      =   find(([pass_folder.isdir]));
indexFiles      =   find(not([pass_folder.isdir]));
nPasses          =   length(indexaDirs);             % number of input files
basePath=pwd;

output_File = 'results_S3A_Dedop_Crete_39.csv';

%fidResults=fopen('./results_CS2_SAR_CMin_FULL2_Svalbard.csv','a+');
% fidResults=fopen('./results_S3_Recheck_6_7_Crete.csv','a+');
%fidResults=fopen('./results_S3_Dedop_Crete_25.csv','a+');
%fidResults=fopen('./results_S3A_Dedop_Crete_36.csv','a+');
%fidResults=fopen('./results_S3B_Dedop_Crete_11_12.csv','a+');
%fidResults=fopen('./results_CR2_20180911_Svalbard_SARIn.csv','a+');
%fidResults=fopen('./results_CR2_20171129_Svalbard_LRM.csv','a+');

fidResults=fopen(['./' output_File],'a+');

for i_pass=3:3
    cd([basePath '\' pass_folder(indexaDirs(i_pass)).name]);
    if(strcmp(mission,'S3'))
        if(i_pass==1)
            fprintf(fidResults,'Cycle;Data; Range error fitting [mm]; Range error minimum value [mm]; Range error aligned ranges [mm]; Datation error fitting [microseconds]; Datation error minimum value[microseconds]; alignment [mm/beam]; noise [mm]; pos; wet[mm]; dry[mm]; iono[mm]; solid_earth[mm]; geocentric_tide[mm]; ocean_loading[mm];geophysical_correction [mm]; alt_rate [m s-1]; zero padding; init beam; final beam; valid beams; L1A File Type; L2 Corrections Type; TRP Internal delay [m]; Track distance [m];\n');
            fclose(fidResults);
        end
    else
        if(i_pass==1)
            fprintf(fidResults,'Data; Range error fitting [mm]; Range error minimum value [mm]; Range error aligned ranges [mm]; Datation error fitting [microseconds]; Datation error minimum value[microseconds]; alignment [mm/beam]; noise [mm]; pos; wet[mm]; wet min pos; dry[mm]; dry min pos; iono[mm]; solid_earth[mm]; geocentric_tide[mm]; ocean_loading[mm];geophysical_correction [mm]; alt_rate [m s-1]; zero padding; init beam; final beam; valid beams; L1A File Type; L2 Corrections Type; TRP Internal delay [m]; Internal delay [m]; Distance[m]; Satellite view;\n');
            fclose(fidResults);
        end
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
    %
    if(strcmp(mission,'S3'))
        clear input_path;
        input_path = dir('./inputs/S3*_SR_1*NT_002.SEN3/');
        mode_L1='NT';
        if(isempty(input_path))
            input_path = dir('./inputs/S3*_SR_1*ST_002.SEN3/');
            mode_L1='ST';
        end
        if(isempty(input_path))
            input_path = dir('./inputs/S3*_SR_1*NT_003.SEN3/');
            mode_L1='NT';
            if(isempty(input_path))
                input_path = dir('./inputs/S3*_SR_1*ST_003.SEN3/');
                mode_L1='ST';
            end
        end
        input_path = strcat(input_path.name, '/');
    end
    %}
    %{
    input_path = dir('./inputs/S3A_SR_1*_002.SEN3');
    if (isempty(input_path))
        input_path = dir('./inputs/S3A_SR_1*_003.SEN3');
    end
    mode_L1='NT';
    input_path = ['/' input_path.name];
    %}
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
    
    if(strcmp(mission,'S3'))
        filesBulk.inputPath       =   strcat('./inputs/',input_path);
    else
        filesBulk.inputPath       =   input_path;
    end
    disp(strcat('Input_path: ',filesBulk.inputPath))
    filesBulk.resultPath      =   output_path;
    disp(strcat('Result_path: ',filesBulk.resultPath))
    mkdir(filesBulk.resultPath);
    mkdir([filesBulk.resultPath 'data/']);
    mkdir([filesBulk.resultPath 'plots/']);
    if(log_flag)
        filesBulk.fid_log = fopen([output_path 'LogError.txt'],'a+');
    end
    filesBulk.inputFiles      =   dir(filesBulk.inputPath);
    filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
    filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
    filesBulk.nFiles          =   length(filesBulk.indexaDirs);             % number of input files
    
    aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
    if(strcmp(mission,'CR2'))
        
        filesBulk.indexFilesL1A=find(~cellfun(@isempty,strfind(aux,'DBL')));
        filesBulk.filterDATAFILES=(~cellfun(@isempty,strfind(aux,'DBL'))) | (~cellfun(@isempty,strfind(aux,'HDR')));
        
        CNF_file=dir('./inputs/cnf_file.m');
        CNF_file=CNF_file.name;
        CHD_file=dir('./inputs/chd_file.m');
        CHD_file=CHD_file.name;
        run(CHD_file);
        run(CNF_file);
        
    elseif(strcmp(mission,'S3'))
        
        filesBulk.indexFilesL1A=find(~cellfun(@isempty,strfind(aux,'nc')));
        filesBulk.filterDATAFILES=(~cellfun(@isempty,strfind(aux,'nc'))) | (~cellfun(@isempty,strfind(aux,'HDR')));
        
        CNF_file=dir('../config/cnf_file.m');
        CNF_file=CNF_file.name;
        CHD_file=dir('../config/chd_file.m');
        CHD_file=CHD_file.name;
        run(CHD_file);
        run(CNF_file);
        
    end
    filesBulk.nFilesL1A=length(filesBulk.indexFilesL1A);
    filesBulk.L1AFiles=filesBulk.inputFiles(filesBulk.indexFilesL1A);
    % for i_fileL1A_input=1:filesBulk.nFilesL1A
    %     filesBulk.L1AFiles(i_fileL1A_input).name    = [filesBulk.L1AFiles(i_fileL1A_input).name '/inputs/measurement_l1a.nc'];
    %     filesBulk.L1AFiles(i_fileL1A_input).isdir   = 0;
    %     listing = dir([filesBulk.inputPath filesBulk.L1AFiles(i_fileL1A_input).name]);
    %     filesBulk.L1AFiles(i_fileL1A_input).bytes   = listing.bytes;
    %
    % end
    
    disp('Total number of L1A to be processed');
    disp(num2str(filesBulk.nFilesL1A))
    
    %---- Display the main processing configuration parameters ------
    
    
    
    %{
    CNF_file=[filesBulk.inputPath filesBulk.inputFiles(~cellfun(@isempty,strfind(aux,'cnf_file'))).name];
    CHD_file=[filesBulk.inputPath filesBulk.inputFiles(~cellfun(@isempty,strfind(aux,'chd_file'))).name];
    %}
    
    
    
    global hamming_window_cnf height_rate_application_cnf FAI_application_cnf CAL2_flag_cnf CAL1p2p_flag_cnf zp_fact_range_cnf
    global apply_stack_mask_cnf avoid_beams_mask_allzeros_cnf avoid_noisy_beams_cnf method_noise_est_cnf param_std_noise_thres_cnf noise_estimation_window
    disp('------------ Main processing configuration parameters -------------')
    disp(strcat('height_rate_application_cnf: ',num2str(height_rate_application_cnf)));
    disp(strcat('FAI_application_cnf: ',num2str(FAI_application_cnf)));
    disp(strcat('CAL2_flag_cnf: ',num2str(CAL2_flag_cnf)));
    disp(strcat('CAL1p2p_flag_cnf: ',num2str(CAL1p2p_flag_cnf)));
    disp(strcat('hamming_window_cnf: ',num2str(hamming_window_cnf)));
    disp(strcat('zp_fact_range_cnf: ',num2str(zp_fact_range_cnf)));
    disp(strcat('apply_stack_mask_cnf: ',num2str(apply_stack_mask_cnf)));
    disp(strcat('avoid_beams_mask_allzeros_cnf: ',num2str(avoid_beams_mask_allzeros_cnf)));
    disp(strcat('avoid_noisy_beams_cnf: ',num2str(avoid_noisy_beams_cnf)));
    if avoid_noisy_beams_cnf
        disp(strcat('method_noise_est_cnf: ',method_noise_est_cnf));
        disp(strcat('param_std_noise_thres_cnf: ',num2str(param_std_noise_thres_cnf)));
        disp(strcat('noise_estimation_window: ',num2str(noise_estimation_window)));
    end
    disp('-------------------------------------------------------------------')
    
    L1B_processing_ACDC_bis(filesBulk,1,0,options);
    Logtime = datestr(now, 'yyyymmddTHHMMSS');
    if(strcmp(mission,'S3'))
        fprintf(filesBulk.fid_log,'%s -> OKKKK!: File %s succesfully processed\n',Logtime);
    else
        fprintf(filesBulk.fid_log,'%s -> OKKKK!: File %s succesfully processed\n',Logtime);%, filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name);
    end
    %{
    for i_fileL1A_input=1:filesBulk.nFilesL1A
        try
            L1B_processing_ACDC_bis(filesBulk,i_fileL1A_input,0,options);
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
    %}
    
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

end


