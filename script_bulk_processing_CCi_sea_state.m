clear;
close all;
% clc;

ttotal=tic;
num_pools=1;
targz_option_active=0;
mission='S3';

%%  
% input_path='/media/alba/DATA/isardSAT/coding/data/CCI_sea_state/L1A_ESA/S3A_SR_1_SRA_A__20200525T053649_20200525T062719_20200526T211407_3029_058_319______MAR_O_ST_004'  
% input_path='/media/alba/DATA/isardSAT/coding/data/CCI_sea_state/L1A_ESA/track_514/S3A_SR_1_SRA_A__20170217T211521_20170217T220550_20170315T124329_3028_014_257______MAR_O_NT_002/'
% input_path='/media/alba/DATA/isardSAT/coding/data/CCI_sea_state/L1A_ESA/track_741/S3A_SR_1_SRA_A__20170225T201722_20170225T210751_20170323T105750_3029_014_370______MAR_O_NT_002/'
% input_path='/media/alba/DATA/isardSAT/coding/data/CCI_sea_state/L1A_ESA/track_513/S3A_SR_1_SRA_A__20170217T202450_20170217T211519_20170315T110727_3029_014_256______MAR_O_NT_002/' 
% input_path='/media/alba/DATA/isardSAT/coding/data/CCI_sea_state/L1A_ESA/track_513/S3A_SR_1_SRA_A__20180329T202459_20180329T211527_20180424T110744_3028_029_256______MAR_O_NT_003/' 
input_path='/media/alba/DATA/isardSAT/coding/data/CCI_sea_state/L1A_ESA/track_513/S3A_SR_1_SRA_A__20191113T202500_20191113T211528_20191209T114804_3027_051_256______MAR_O_NT_003/' 

output_path='/media/alba/DATA/isardSAT/coding/output/CCI_sea_state/ACDC/';
dirs=split(input_path,'/');
mkdir([[output_path '/'] char(dirs(end-1))]);
output_path=fullfile(output_path,char(dirs(end-1)))

GPP_bulk_processing_paralelization(mission,input_path,output_path,num_pools,targz_option_active)

time=toc(ttotal);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time (all regions): ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);
