global generate_kml
%% ------------------------ GENERAL COMMON SETTINGS -----------------------
code_folder_full_path='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/Matlab_DeDop_code/';
figures_visible=0;
flag_outliers_removal=1;
type_outliers_removal='hampel';
active_validation_tracks=1;
active_comparison_tracks=0;
num_pools=1;
sh_name_nc='ssh'; %name of the variable in netcdf refering to the surface height

%-------------------------- INCLUDING SEARCH PATH CODE --------------------
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


%% ------------------------ CENTRAL PACIFIC ----------------------------------------------------------------------------
% A la cryosat2
disp('--------------------------------------------------------------------')
disp('-------------------- CENTRAL PACIFIC ----------------------------------')
disp('--------------------------------------------------------------------')
input_path_L2_ISD='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/results/2013_L1B_ISD_CR2/central_pacific/data/'
input_path_L2_ESA=''
path_comparison_results='C:/Users/eduard.makhoul/isardSAT/conferences_meetings/OSTST2016/processing/data/results/L1B_ACDC_CR2/central_pacific/validation_new/'
filename_mask_KML='';
try
    generate_kml=1;
    validation_L1B_ACDC_bulk_region_parallel(input_path_L2_ISD,path_comparison_results,...
                                            'num_pools',num_pools,...
                                            'figures_visible',figures_visible,...
                                            'flag_outliers_removal',flag_outliers_removal,...
                                            'type_outliers_removal',type_outliers_removal,...
                                            'active_comparison_tracks',active_comparison_tracks,...
                                            'active_validation_tracks',active_validation_tracks,...
                                            'filename_mask_KML',filename_mask_KML,...
                                            'input_path_L2_ESA',input_path_L2_ESA,...
                                            'sh_name_nc',sh_name_nc);
catch
    disp('Some errors occurred to complete the whole validation')
end
