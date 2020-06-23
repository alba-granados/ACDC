    clear;
close all;
% clc;

ttotal=tic;
num_pools=1;
targz_option_active=0;
mission='S3';

%%  
input_path='/media/alba/DATA/isardSAT/coding/data/CCI_sea_state/L1A_ESA/S3A_SR_1_SRA_A__20200525T053649_20200525T062719_20200526T211407_3029_058_319______MAR_O_ST_004'  
output_path='/media/alba/DATA/isardSAT/coding/output/CCI_sea_state/ACDC'
GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)

time=toc(ttotal);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time (all regions): ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);
