clear;
close all;
clc;

ttotal=tic;
num_pools=1;
targz_option_active=0;
mission='CR2';

%% ROI_3 central pacific
input_path='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/FBR_ESA/2013/central_pacific/subset/'
output_path='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/results/2013_L1B_ISD_CR2/central_pacific/subset/'
GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)

% %% ROI_1 west pacific
% input_path='/data/DATA/SCOOP/disk_copy/west_pacific/'
% output_path='/data/DATA/SCOOP/results/L1B_ISD/west_pacific/'
% GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)
% 
% %% ROI_2 Central Pacific
% input_path='/data/DATA/SCOOP/disk_copy/central_pacific/'
% output_path='/data/DATA/SCOOP/results/L1B_ISD/central_pacific/'
% GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)
% 
% %% ROI_5 North Sea
% input_path='/data/DATA/SCOOP/disk_copy/north_sea/'
% output_path='/data/DATA/SCOOP/results/L1B_ISD/north_sea/'
% GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)

% %% ROI_6 Agulhas
% input_path='/data/DATA/SCOOP/copy_disk_ESA_2013/agulhas/'
% output_path='/data/DATA/SCOOP/results/L1B_ISD/agulhas/2013/'
% GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)
% 
% 
% %% ROI_8 North Indian Coast
% input_path='/data/DATA/SCOOP/copy_disk_ESA_2013/north_indian_coast/'
% output_path='/data/DATA/SCOOP/results/L1B_ISD/north_indian_coast/2013/'
% GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)
% 
% %% ROI_9 Indonesia
% input_path='/data/DATA/SCOOP/copy_disk_ESA_2013/indonesia/'
% output_path='/data/DATA/SCOOP/results/L1B_ISD/indonesia/2013/'
% GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)
% 
% %% ROI_10 Harvest
% input_path='/data/DATA/SCOOP/copy_disk_ESA_2013/harvest/'
% output_path='/data/DATA/SCOOP/results/L1B_ISD/harvest/2013/'
% GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)

time=toc(ttotal);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time (all regions): ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);