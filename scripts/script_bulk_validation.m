%% ------------------------ Setting folders search ------------------------
%code_folder_full_path='/home/emak/feina/projects/SCOOP/processing/Matlab_DeDop_code/';
code_folder_full_path='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/Matlab_DeDop_code/';
cd(code_folder_full_path);
FolderInfo=dir(code_folder_full_path);
FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
aux=struct2cell(FolderInfo); %into a cell array where first row is 
folders=(aux(1,[FolderInfo.isdir]))'; %name is the first row
clear aux;
folders=strcat(code_folder_full_path,folders(~strcmp(folders,['.'])&~strcmp(folders,['..'])&~strcmp(folders,['.svn'])&~strcmp(folders,['inputs'])));
for i_folder=1:length(folders)
    addpath(genpath(char(folders(i_folder))));
end
set_default_plot

num_pools=1;

%% ------------------------ AMAZON ----------------------------------------------------------------------------
disp('--------------------------------------------------------------------')
disp('-------------------------- AMAZON ----------------------------------')
disp('--------------------------------------------------------------------')
% input_path_L1B_ISD='/data/DATA/SCOOP/results/L1B_ISD/west_pacific/data/'
% input_path_L1B_ESA='/data/DATA/SCOOP/ESA/L1B/west_pacific/'
% input_path_L2_ESA='/data/DATA/SCOOP/ESA/L2/west_pacific/'
% path_comparison_results='/data/DATA/SCOOP/results/validation_L1B_S3/west_pacific/'
input_path_L1B_ISD='C:/Users/eduard.makhoul/isardSAT/projects/SHAPE/data/results/L1B_ISD_S3/Amazon/data/subset/'
input_path_L1B_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SHAPE/data/input/L1B_ESA/Amazon/'
input_path_L2_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SHAPE/data/input/L2_ESA/Amazon/'
path_comparison_results='C:/Users/eduard.makhoul/isardSAT/projects/SHAPE/data/results/validation/Amazon/subset/'
validation_bulk_region_EM_parallel(input_path_L1B_ISD,input_path_L1B_ESA,input_path_L2_ESA,path_comparison_results,'num_pools',num_pools);


% %% ------------------------ WEST PACIFIC ----------------------------------------------------------------------------
% disp('--------------------------------------------------------------------')
% disp('-------------------- WEST PACIFIC ----------------------------------')
% disp('--------------------------------------------------------------------')
% % input_path_L1B_ISD='/data/DATA/SCOOP/results/L1B_ISD/west_pacific/data/'
% % input_path_L1B_ESA='/data/DATA/SCOOP/ESA/L1B/west_pacific/'
% % input_path_L2_ESA='/data/DATA/SCOOP/ESA/L2/west_pacific/'
% % path_comparison_results='/data/DATA/SCOOP/results/validation_L1B_S3/west_pacific/'
% input_path_L1B_ISD='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L2_processor_validation/input_data/L1B/ISD_nc/west_pacific_CR2/'
% input_path_L1B_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L1B_ESA/west_pacific/'
% input_path_L2_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L2_ESA/west_pacific/'
% path_comparison_results='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L2_processor_validation/output_results/west_pacific_CR2/plots/others/'
% validation_bulk_region_EM_parallel(input_path_L1B_ISD,input_path_L1B_ESA,input_path_L2_ESA,path_comparison_results,'num_pools',num_pools);

% %% ------------------------ CENTRAL PACIFIC ----------------------------------------------------------------------------
% disp('--------------------------------------------------------------------')
% disp('----------------- CENTRAL PACIFIC ----------------------------------')
% disp('--------------------------------------------------------------------')
% % input_path_L1B_ISD='/data/DATA/SCOOP/results/L1B_ISD/central_pacific/data/'
% % input_path_L1B_ESA='/data/DATA/SCOOP/ESA/L1B/central_pacific/'
% % input_path_L2_ESA='/data/DATA/SCOOP/ESA/L2/central_pacific/'
% % path_comparison_results='/data/DATA/SCOOP/results/validation_L1B_S3/central_pacific/'
% input_path_L1B_ISD='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/results/L1B_ISD/central_pacific_CR2/'
% input_path_L1B_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L1B_ESA/central_pacific/'
% input_path_L2_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L2_ESA/central_pacific/'
% path_comparison_results='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L1B_validation/A_la_CryoSat_2/central_pacific/'
% validation_bulk_region_EM_parallel(input_path_L1B_ISD,input_path_L1B_ESA,input_path_L2_ESA,path_comparison_results,'num_pools',num_pools);

% %% ------------------------ EAST PACIFIC ----------------------------------------------------------------------------
% disp('--------------------------------------------------------------------')
% disp('-------------------- EAST PACIFIC ----------------------------------')
% disp('--------------------------------------------------------------------')
% % input_path_L1B_ISD='/data/DATA/SCOOP/results/L1B_ISD/east_pacific/data/'
% % input_path_L1B_ESA='/data/DATA/SCOOP/ESA/L1B/east_pacific/'
% % input_path_L2_ESA='/data/DATA/SCOOP/ESA/L2/east_pacific/'
% % path_comparison_results='/data/DATA/SCOOP/results/validation_L1B_S3/east_pacific/'
% input_path_L1B_ISD='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/results/L1B_ISD_S3/east_pacific/subset/'
% input_path_L1B_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L1B_ESA/east_pacific/'
% input_path_L2_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L2_ESA/central_pacific/'
% path_comparison_results='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L1B_validation/A_la_CryoSat_2/central_pacific/'
% validation_bulk_region_EM_parallel(input_path_L1B_ISD,input_path_L1B_ESA,input_path_L2_ESA,path_comparison_results,'num_pools',num_pools);
% 
% %% ------------------------ HARVEST ----------------------------------------------------------------------------
% disp('--------------------------------------------------------------------')
% disp('------------------------- HARVEST ----------------------------------')
% disp('--------------------------------------------------------------------')
% % input_path_L1B_ISD='/data/DATA/SCOOP/results/L1B_ISD/central_pacific/data/'
% % input_path_L1B_ESA='/data/DATA/SCOOP/ESA/L1B/central_pacific/'
% % input_path_L2_ESA='/data/DATA/SCOOP/ESA/L2/central_pacific/'
% % path_comparison_results='/data/DATA/SCOOP/results/validation_L1B_S3/central_pacific/'
% input_path_L1B_ISD='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/results/L1B_ISD_S3/harvest/'
% input_path_L1B_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L1B_ESA/harvest/'
% input_path_L2_ESA='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L2_ESA/harvest/'
% path_comparison_results='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L1B_validation/A_la_Sentinel_3/harvest/'
% validation_bulk_region_EM_parallel(input_path_L1B_ISD,input_path_L1B_ESA,input_path_L2_ESA,path_comparison_results,'num_pools',num_pools);
