% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT Ltd. 
% ---------------------------------------------------------
% Objective: Retreive input data from the inputs folder. Read auxiliaries
% ---------------------------------------------------------
% 
% Author:   Albert Garcia-Mondejar / isardSAT
%           Roger Escola Jane      / isardSAT
%
% Version  record
% 1.0 2015/01/01 First Version
% 2.0 2016/03/22 Removed the reading routine for the science. Added meteo reading for CR2 
% 2.1 2016/05/20 Added case when there are no records inside the mask
% 2.2 2016/06/01 S3 case corrections
% 2.3 2016/06/04 global product_coord created to be used when plotting full
% track of the product or selected records within the mask
% 2.4 2016/11/25 Included plot_track and write kml with GUI options. Also
% added exit_gpp as output when no reco
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [steps,meteo,files,first_burst,final_burst] = read_inputs(files)

global mode level mission product_coord

steps = [1,1,1,0];
% steps for the processor
% steps(1) -> L0    to L1A
% steps(2) -> L1A   to L1B-S
% steps(3) -> L1B-S to L1B
% steps(4) -> L1B   to L2

i_ISP_file  = 1;

meteo = [];
files.characterization=[];
for iFile   =   1:files.nFiles
    iFilename           = files.inputFiles(files.indexFiles(iFile)).name;
    if(strfind(iFilename,'chd_file'))
        files.filename_CHD    = iFilename;
    elseif(strfind(iFilename,'cst_file'))
        files.filename_CST    = iFilename;
    elseif(strfind(iFilename,'cnf_file'))
        files.filename_CNF    = iFilename;
    elseif(strfind(iFilename,'rep_file'))
        files.filename_REP    = iFilename;
    elseif(strfind(iFilename,'spw_file'))    
        files.filename_SPW    = iFilename;
    elseif(strfind(iFilename,'FAT'))    
        files.filename_FAT    = iFilename;
    elseif(strfind(iFilename,'FSV'))    
        files.filename_OSV    = iFilename;
    elseif(strfind(iFilename,'ISP'))    
        files.filename_ISP(i_ISP_file,:)  = iFilename;
%        size(i_ISP_file)            = files.inputFiles(files.indexFiles(iFile)).bytes;
        i_ISP_file      = i_ISP_file+1;
        mission    = 'S6_';
        level      = 'L0';
    elseif(strfind(iFilename,'measurement_l1a'))
        files.filename_L1A    = [files.inputPath iFilename];
        hdr     = dir([files.inputPath '*.xml']);
        files.filename_xmldump = [files.inputPath hdr.name];
        files.xfdumanifest = xmlread(files.filename_xmldump);
        i_ISP_file      = i_ISP_file+1;
        mission    = 'S3_';
        level      = 'L1A';
		mode       = 'SAR';
    elseif(strfind(iFilename,'SIN_FR'))    
        if (strfind(iFilename,'.DBL'))
            files.filename_L1A    = [files.inputPath iFilename];
            size_L1A        = files.inputFiles(files.indexFiles(iFile)).bytes;
            mission    = 'CR2';
            level      = 'L1A';
            mode       = 'SIN';            
            files.sph = readSPH(files.filename_L1A);
            fbr_dsd = files.sph.dsds(1);
            files.ds_offset = fbr_dsd.ds_offset;
            files.dsr_size = fbr_dsd.dsr_size;
        end
    elseif(strfind(iFilename,'SAR_FR'))   
        if (strfind(iFilename,'.DBL'))
            files.filename_L1A    = [files.inputPath iFilename];
            size_L1A        = files.inputFiles(files.indexFiles(iFile)).bytes;
            mission    = 'CR2';
            level      = 'L1A';
            mode       = 'SAR';                       
            files.sph = readSPH(files.filename_L1A);
            fbr_dsd = files.sph.dsds(1);
            files.ds_offset = fbr_dsd.ds_offset;
            files.dsr_size = fbr_dsd.dsr_size;
        end
        %
    elseif(strfind(iFilename,'*.kml'))
        files.filename_kml    = [files.inputPath iFilename];
        %}
    elseif(strfind(iFilename,'characterization'))
        %Nc file containing the calibration files
        files.characterization = [files.inputPath iFilename];
    end
    if(strfind(iFilename,'.kml'))
        files.filename_kml    = [files.inputPath iFilename];
    end
    
end

iConSize = size(files.configPath);
iConSize = iConSize(1);
%mission = 'S3';
%if (strcmp(mission,'S3_')==1||strcmp(mission,'S3')==1)
    for iCon = 3:iConSize
        iConfigfile         = files.configPath(iCon);
        if(strfind(iConfigfile.name,'chd_file'))
            files.filename_CHD    = iConfigfile;
        elseif(strfind(iConfigfile.name,'cst_file'))
            files.filename_CST    = iConfigfile;
        elseif(strfind(iConfigfile.name,'cnf_file'))
            files.filename_CNF    = iConfigfile;
        elseif(strfind(iConfigfile.name,'spw_file'))
            files.filename_SPW    = iConfigfile;
        elseif(strfind(iConfigfile.name,'.kml'))
            files.filename_kml    = iConfigfile;
        end
    end
%end
%% 1.  UNPACKING & DECODING INPUTS

if (strcmp(mission,'S3_')==1||strcmp(mission,'S3')==1)  % why? 
    run([files.configDir files.filename_CST.name]);
    run([files.configDir files.filename_CHD.name]);
    run([files.configDir files.filename_CNF.name]);
    run([files.configDir files.filename_SPW.name]);
else
    run([files.inputPath files.filename_CST]);
    run([files.inputPath files.filename_CHD]);
    run([files.inputPath files.filename_CNF]);
    run([files.inputPath files.filename_SPW]); 
end
    global mask_flag
    global N_bursts_cycle_chd N_bursts_cycle_sar_chd N_bursts_cycle_sarin_chd N_ku_pulses_burst_chd
    global burst_phase_array_cor_cal1_sar_rep_1 burst_phase_array_cor_cal1_sar_rep_2 burst_power_array_cor_cal1_sar_rep_1 burst_power_array_cor_cal1_sar_rep_2 
    global wfm_cal2_science_sar_rep_1 wfm_cal2_science_sar_rep_2
    global N_samples N_samples_sar_chd N_samples_sin_chd zp_fact_range_cnf
    global bri_chd brf_chd bri_sar_chd brf_sar_chd bri_sin_chd brf_sin_chd
    global window_rg_cnf window_range_SAR window_range_SARin window_range
    global gain_scale_win_rg_SAR gain_scale_win_rg_SARin gain_scale_win_rg
    global onboard_proc_sar_chd
    global bw_ku_chd pulse_length_chd
  % alba: added antenna_beamwidth from chd_file to be used in nf_p.alphax (gen_nonfit_params_EM.m) to compute Bkl (stack_gen.m) 
    global antenna_beamwidth_ku_chd antenna_gain_ku_chd
    global antenna_beamwidth_alt_ku_chd antenna_beamwidth_act_ku_chd
    
    switch mode
        case 'SAR'
            N_samples = N_samples_sar_chd;    
            N_bursts_cycle_chd = N_bursts_cycle_sar_chd;
            bri_chd=bri_sar_chd;
            brf_chd=brf_sar_chd;
            window_range=window_range_SAR;
            gain_scale_win_rg=gain_scale_win_rg_SAR;

        case 'SIN'                    
            N_samples = N_samples_sin_chd;
            N_bursts_cycle_chd = N_bursts_cycle_sarin_chd;
            bri_chd=bri_sin_chd;
            brf_chd=brf_sin_chd;
            window_range=window_range_SARin;
            gain_scale_win_rg=gain_scale_win_rg_SARin;                    
    end
            
   
    if ~isempty(files.characterization)                
        %read the calibration files
        switch mode
            case 'SAR'
                burst_phase_array_cor_cal1_sar_rep_1 = ncread(files.characterization,'cal1_p2p_phase_sar').'; % in radians
                burst_power_array_cor_cal1_sar_rep_1 = 20*log10(ncread(files.characterization,'cal1_p2p_amplitude_sar').'); % in dB
                wfm_cal2_science_sar_rep_1 = ncread(files.characterization,'cal2_lpf_sar').';
            case 'SIN'
                burst_phase_array_cor_cal1_sar_rep_1 = ncread(files.characterization,'cal1_p2p_phase_sarin_rx1').'; % in radians
                burst_phase_array_cor_cal1_sar_rep_2 = ncread(files.characterization,'cal1_p2p_phase_sarin_rx2').'; % in radians
                burst_power_array_cor_cal1_sar_rep_1 = 20*log10(ncread(files.characterization,'cal1_p2p_amplitude_sarin_rx1').'); % in dB
                burst_power_array_cor_cal1_sar_rep_2 = 20*log10(ncread(files.characterization,'cal1_p2p_amplitude_sarin_rx2').'); % in dB
                wfm_cal2_science_sar_rep_1 = ncread(files.characterization,'cal2_lpf_sarin_rx1').';
                wfm_cal2_science_sar_rep_2 = ncread(files.characterization,'cal2_lpf_sarin_rx2').';
        end

    else
        switch mode
            case 'SAR'
                burst_phase_array_cor_cal1_sar_rep_1 =zeros(1,N_ku_pulses_burst_chd); % in radians
                burst_power_array_cor_cal1_sar_rep_1 = zeros(1,N_ku_pulses_burst_chd); % in dB
                wfm_cal2_science_sar_rep_1 = ones (1,N_samples);
            case 'SIN'
                burst_phase_array_cor_cal1_sar_rep_1 = zeros(1,N_ku_pulses_burst_chd); % in radians
                burst_phase_array_cor_cal1_sar_rep_2 = zeros(1,N_ku_pulses_burst_chd); % in radians
                burst_power_array_cor_cal1_sar_rep_1 = zeros(1,N_ku_pulses_burst_chd); % in dB
                burst_power_array_cor_cal1_sar_rep_2 = zeros(1,N_ku_pulses_burst_chd); % in dB
                wfm_cal2_science_sar_rep_1 = ones (1,N_samples);
                wfm_cal2_science_sar_rep_2 = ones (1,N_samples);
        end
    end
    
    switch level
        case 'L0'
            steps(1)    = 1;
            switch mission
                case 'CR2'
                    switch mode
                        case 'SAR'
                        case 'SIN'
                            steps(1)    = 0;
                            %                             [input] = readFBR_CryoSat_SIN(files.filename_L1A,size_L1A);
                    end
                case 'S3_'
                    [netCDF_L0] = readanyNETCDF_V1(files.filename_ISP);
                    %                             [input]     = adapt_S3L0(netCDF_L0);
                    
                    
                    %                             [meteo]     = read_meteo_rinex(files.filename_meteo);
                    
                case 'S6_'
                    %                             [input]= read_ISP(filename_ISP(i_ISP,:), size_ISP(i_ISP));
                    
            end
        case 'L1A'
            steps(1)    = 0;
            switch mission
                case 'S6_'
                case {'S3_','S3'}
                    
                    if (mask_flag)
%                         product_mask  = kml2lla(strcat(files.configDir,files.filename_kml.name));
                        product_mask  = kml2lla(files.filename_kml);
                        [lat,lon,alt,range_delay] = readnetCDF_Sentinel3_lat_lon(files.filename_L1A);
                        alt_surf= alt - range_delay;
                        index_inside = find(inpolygon(lon,lat,product_mask.coord(:,1),product_mask.coord(:,2))); % records inside the given mask
                        if(isempty(index_inside))
							first_burst=0;
							final_burst=0;
                        else
                            ncid = netcdf.open(files.filename_L1A,'NC_NOWRITE'); %open file
                            [ndims,~,~,~] = netcdf.inq(ncid); %get global attributes
                            for i_dim=0:(ndims-1)
                                [dimname, dimlen] = netcdf.inqDim(ncid,i_dim);
                                if(strcmp(dimname,'time_l1a_echo_sar_ku'))
                                    break;
                                end
                            end
                            netcdf.close(ncid);
                            first_burst = max(1,(index_inside(1))-N_ku_pulses_burst_chd*N_bursts_cycle_chd);
                            final_burst = ((max(index_inside)+N_ku_pulses_burst_chd*N_bursts_cycle_chd));
                            if(isempty(strfind(files.filename_L1A,'measurement_l1a_reduced')))
                                [files] = create_NetCDF_L1A_S3(files,final_burst-first_burst);
                                resize_NetCDF_L1A_S3(files.filename_L1A,first_burst, final_burst-first_burst, files); % alba: added , files
%                                 movefile('./results/data/measurement_l1a_reduced.nc','./inputs/');
                                movefile(strcat(strcat(files.resultPath,'data/'),'measurement_l1a_reduced.nc'),files.inputPath) % alba
                                %movefile(files.filename_L1A,'./results/data/');
                                final_burst=final_burst-first_burst;
                                first_burst=1;
%                                 files.filename_L1A='./inputs/measurement_l1a_reduced.nc';
                                files.filename_L1A=strcat(files.inputPath,'measurement_l1a_reduced.nc');
                                
                            end
						end
                    else
                        
                        if(files.options.plotting_flag(3) || files.options.writting_flag(4))
                           [lat,lon,alt,range_delay] = readnetCDF_Sentinel3_lat_lon(files.filename_L1A);
                            alt_surf= alt - range_delay; 
                        end
                        ncid = netcdf.open(files.filename_L1A,'NC_NOWRITE'); %open file
                        [ndims,~,~,~] = netcdf.inq(ncid); %get global attributes
                        for i_dim=0:(ndims-1)
                            [dimname, dimlen] = netcdf.inqDim(ncid,i_dim);
                            if(strcmp(dimname,'time_l1a_echo_sar_ku'))
                                break;
                            end
                        end
                        netcdf.close(ncid);
                        first_burst = 1;
                        final_burst = dimlen;
                    end
                    
                case 'CR2'                    
                    if (mask_flag)
                        product_mask  = kml2lla(files.filename_kml);
                        [lat,lon,alt,win_delay,~] = readFBR_CryoSat_SARSIN_lat_lon(files.filename_L1A);
                        alt_surf= alt - win_delay * c_cst/2;
                        index_inside = find(inpolygon(lon,lat,product_mask.coord(:,1),product_mask.coord(:,2))); % records inside the given mask
						
						if(isempty(index_inside))
							first_burst=0;
							final_burst=0;
						else
							first_burst = max(1,(index_inside(1))-N_ku_pulses_burst_chd*N_bursts_cycle_chd);
							final_burst = min((max(index_inside)+N_ku_pulses_burst_chd*N_bursts_cycle_chd),files.sph.dsds(1).num_dsr*20);
							product_coord = [lat(first_burst:final_burst); lon(first_burst:final_burst); alt(first_burst:final_burst); alt_surf(first_burst:final_burst)];

						end
                        %                                 N_bursts = final_burst -first_burst +1;
                    else
                        if(files.options.plotting_flag(3) || files.options.writting_flag(4))
                            [lat,lon,alt,win_delay,~] = readFBR_CryoSat_SARSIN_lat_lon(files.filename_L1A);
                            alt_surf= alt - win_delay * c_cst/2;
                        end
                        first_burst = 1;
                        final_burst = (files.sph.dsds(1).num_dsr-1)*20;
                        
                    end
                    %onboard_proc_sar_chd=onboard_proc_sar_chd+10*log10(N_samples)+10*log10(N_ku_pulses_burst_chd);
                    %onboard_proc_sar_chd=onboard_proc_sar_chd+10*log10(N_samples*zp_fact_range_cnf);
                    %second option using the TBP
                    %onboard_proc_sar_chd =onboard_proc_sar_chd+10*log10(bw_ku_chd*pulse_length_chd); % Processing gain only on pulse compression: TBP
                    %                             [meteo]     = read_CryoSat_GeoCorr_rinex(files.filename_L1A,first_record, num_records);
            end
            if(files.options.plotting_flag(3))
                [product_mask] = plot_track(lat,lon,alt_surf,files.options,product_mask);
                if (mask_flag)
                    index_inside = find(inpolygon(lon,lat,product_mask.coord(:,1),product_mask.coord(:,2))); % records inside the given mask
                    if(isempty(index_inside))
                        first_burst=0;
                        final_burst=0;
                    else
                        first_burst = max(1,(index_inside(1))-N_ku_pulses_burst_chd*N_bursts_cycle_chd);
                        final_burst = ((max(index_inside)+N_ku_pulses_burst_chd*N_bursts_cycle_chd));

                    end                    
                end
                
            end
            if(files.options.writting_flag(4))
                disp('l1a2kml...')
                lla2kml(files.filename_L1A,lat,lon,alt_surf,'.');
                disp('Done l1a2kml.')
                movefile([files.filename_L1A '.kml'],files.resultPath);
            end
    end

end
