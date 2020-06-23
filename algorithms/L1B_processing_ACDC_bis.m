% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT Ltd.
% --------------------------------------------------------
%
% Ground Processor Prototype for altimetry missions:
% CryoSat-2 SARIn
% CryoSat-2 SAR
% Sentinel 6 RAW
% Sentinel 6 RMC
%
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
% Authors:  Albert Garcia-Mondejar  / isardSAT
%           Eduard Mackoul          / isardSAT
%           Roger Escola Jane       / isardSAT
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v1.0 2016/04/xx   First version based on S6 GPP
% v2.0 2016/04/xx   Changed reading routines for each record
% v2.1 2016/04/xx   Changed for bulk processing. Files are opened one time within read_L1A_record
% v2.2 2016/05/12   Changed FBR for L1A. Opening file moved inside read_L1A_record fuction
% v2.3 2016/05/17   Include the possibility to activate or deactivate the alignment of the L1B waveforms w.r.t first one
%					Stacking of the remaining last surfaces
% v2.4 2016/05/18 	Add writting L1BS
% v2.5 2016/05/20 	Inlucded return to skip processing when no records are found within the mask
% v2.6 2016/05/20   Included extend_window_cnf option
% v2.7 2017/02/04   Added log to write progress each 1%

function L1B_processing_ACDC_bis(filesBulk,i_fileL1A_input,targz_option_active, options)
%   writting_L1B = 0;
% 	writting_L1BS = 0;
exit_gpp = 0;

global extend_window_cnf include_wfms_aligned N_ku_pulses_burst_chd
global N_bursts_cycle_chd mission
global optional_ext_file_flag file_ext_string
optional_ext_file_flag=0; % alba: "Flag indicating whether to include or not a given extension to the output product name of file"
file_ext_string='_height_rate';

%EM: added 03.10.2016
global ACDC_application_cnf cnf_p_ACDC figures_visible_cnf
global bw_ku_chd zp_fact_range_cnf N_samples c_cst
if(strcmp(mission,'S3'))
    if (figures_visible_cnf || options.GUI_flag)
        set(0, 'DefaultFigureVisible', 'on');
    else
        set(0, 'DefaultFigureVisible', 'off');
    end
end
t6 = tic;
%% Define the files structure for the specific FBR to be processed
% including all other files (cnf,cst,chd,....) not DBL nor HDR related to the specific
% folder (in order no to change the functions)
files.inputPath =filesBulk.inputPath;
files.configPath=dir('../config/'); % why? I added this folder with cnf, cst, chd files
files.configDir =('../config/');
files.resultPath =filesBulk.resultPath;
files.inputFiles =filesBulk.inputFiles(~(filesBulk.filterDATAFILES));
files.options = options;
if targz_option_active
    %untar the file
    
    untar([files.inputPath filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name],files.inputPath);
    filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name=strrep(filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name,'TGZ','DBL');
end
%add the current L1A to be processed
files.inputFiles(end+1) =filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input));
files.indexaDirs      =   find(([files.inputFiles.isdir]));
files.indexFiles      =   find(not([files.inputFiles.isdir]));
files.nFiles          =   length(files.indexFiles);

filename_L1A=filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name;
disp(strcat('Processing L1A: ',filename_L1A))
disp('....')

%progressbar('Burst processing','Surfaces computation','Burst focussing','Stack processing');
i_surf              = 1;
i_burst_focussed    = 1;
i_surf_stacked      = 1;
% previous_gap_flag   = 0; % indicates that the burst has a gap unable to be fullfilled
%% Read inputs
[steps,meteo,files,first_burst,final_burst] = read_inputs(files);

N_bursts    = final_burst-first_burst+1;
if(N_bursts < N_bursts_cycle_chd*N_ku_pulses_burst_chd)
    disp(['Not enough records inside mask to create a complete stack for the file ' files.filename_L1A]);
    fprintf(filesBulk.fid_log,'%s\n' ,['Not enough records inside mask to create a complete stack for the file ' files.filename_L1A]);
    return;
end
original_burst_index    =   first_burst:final_burst;
burst_indexes_in_memory = 1: N_bursts;
N_surfs_loc_estimated = floor(N_bursts/N_bursts_cycle_chd*2.0);
disp(strcat('# Estimated surfaces',{' '},num2str(N_surfs_loc_estimated)))
if(options.writting_flag(2))
    disp('creating NetCDF_L1B_S3 in output folder...')
    [files] = create_NetCDF_L1B_S3(files,N_surfs_loc_estimated);
end
if(options.writting_flag(1))
    disp('creating NetCDF_L1BS_S3 in output folder...')
    [files] = create_NetCDF_L1BS_S3(files,N_surfs_loc_estimated);
end
[L1A]           = create_L1A_struct;
[L1A_buffer]    = create_L1A_struct;
[L1BS_buffer]   = create_L1BS_struct;
L1B             = [];

% i_burst=1;
load(strcat(files.resultPath,'/data/up_to_i_surf96.mat'))
[steps,meteo,files,first_burst,final_burst] = read_inputs(files);

N_bursts_original=N_bursts;
fprintf('N_bursts=%d\n', N_bursts_original)
while(i_burst<=N_bursts)
    fprintf('\ni_burst=%d\n', i_burst)
    if((floor(i_burst/N_bursts*100))>(floor((i_burst-1)/N_bursts*100)))
        %             fprintf(filesBulk.fid_log,'%s\n' ,['L1A record ' num2str(i_burst) ' of  ' num2str(N_bursts) ': ' num2str(floor(i_burst/N_bursts*100)) '%']);
    end
    %progressbar(i_burst/N_bursts,[],[],[]);
    %% BURST Processor
    if (steps(1))
        fprintf('L0 to L1A processing (step(1)=1)...\n')
        [ISP]  = read_ISP_record(files,first_burst,original_burst_index(i_burst),N_bursts); %To be build
        [L1A]  = pre_datation    (ISP);
        [L1A]  = pre_win_delay   (L1A,ISP);
        [L1A]  = final_burst_datation   (L1A,ISP);
        [L1A]  = onboard_reversion   (L1A,ISP);
        [L1A]  = instrument_gain_corr   (ISP,L1A);
        [L1A]  = waveforms_correction   (ISP,L1A);
    else
        fprintf('Reading L1A input files (skip L0 to L1A processing)...\n')
        [L1A,files]  = read_L1A_record(files,L1A,original_burst_index(i_burst),i_burst,N_bursts);
        if(i_burst > 1)
            %if there is are missing burst records in the L1A file, jumps
            %in time, they are created to provide continuity to the stacks.
            [L1A,original_burst_index,N_bursts] = check_continuity (L1A_buffer(end),L1A,original_burst_index,i_burst,N_bursts);
        end        
        if(L1A.confi_block_degraded==1&&L1A.days==0) % handle degraded bursts
            [L1A,end_of_file]  = fill_empty_burst(files,L1A_buffer(end),L1A,original_burst_index(i_burst),N_bursts_original,i_burst);
            if(end_of_file==1)
                N_bursts=i_burst-1;
                break;
            end
        end        
    end
    L1A_buffer(i_burst) = L1A;
    
    %checking purposes
    %         burst_average_wfm(i_burst,:)=nanmean((abs(fftshift(fft([L1A_buffer(i_burst).wfm_cal_gain_corrected],[],2),2))).^2,1);
    %         win_delay_sar_ku(i_burst)=L1A_buffer(i_burst).win_delay_sar_ku;
    %         cai_sar_isp(i_burst)=L1A_buffer(i_burst).cai_sar_isp;
    %         fai_sar_isp(i_burst)=L1A_buffer(i_burst).fai_sar_isp;
    %         alt_rate_sar_sat(i_burst)=L1A_buffer(i_burst).alt_rate_sar_sat;
    
    
    
    if(steps(2))
        fprintf('Step(2) L1A to L1BS processing...\n')
        fprintf('Computing surface locations...\n')
        [L1BS_buffer, out_surf] = surface_locations (L1A_buffer,L1BS_buffer, i_burst, i_surf);
        if(i_surf < out_surf)
            %progressbar([],i_surf/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd,[],[]);
            i_surf = out_surf;
            N_total_surf_loc = i_surf-1;
            new_surface=1;
            if(files.options.plotting_flag(3))
                plotm(L1BS_buffer(i_surf-1).lat_surf,L1BS_buffer(i_surf-1).lon_surf,'r.','Parent',options.axes);
                plot(L1BS_buffer(i_surf-1).alt_surf, L1BS_buffer(i_surf-1).lat_surf,'r.','Parent',options.wd_axes);
                pause(0.01);
            end
        end
        
%         i_surf = 96;
        fprintf('i_surf=%d\n', i_surf)
%         load(strcat(files.resultPath,'/data/up_to_i_surf96.mat'))
        if(i_surf > 64) %  > SAR Ku pulses in burst (number of surface location that are “observed” by the satellite at the current satellite burst position)
            if(i_burst_focussed==1)% handle final bursts
                burst_margin = i_burst;
            end
            fprintf('beam angles computation and azimuth processing...\n')
            %progressbar([],[],i_burst_focussed/N_bursts,[]);
            [L1A_buffer(i_burst_focussed)]      = beam_angles (L1A_buffer(i_burst_focussed),   L1BS_buffer, N_total_surf_loc,i_surf_stacked);
            [L1A_buffer(i_burst_focussed)]      = azimuth_processing(L1A_buffer(i_burst_focussed));
            i_burst_focussed = i_burst_focussed+1;
            
            fprintf('L1A_buffer(i_burst_focussed-1).surf_loc_index(1)=%d\n', L1A_buffer(i_burst_focussed-1).surf_loc_index(1))
            fprintf('i_burst_focussed=%d\n', i_burst_focussed)
            fprintf('i_surf_stacked=%d\n', i_surf_stacked)
            if(L1A_buffer(i_burst_focussed-1).surf_loc_index(1)~=i_surf_stacked)
%                 save(strcat(files.resultPath,'/data/up_to_i_surf96_2.mat'));
                fprintf('stacking, geometry corrections, range transformation, stack masking...\n')
                %             if(i_burst_focussed > (N_bursts_cycle_chd*N_ku_pulses_burst_chd))
                %progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
                    [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
                [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
                [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
                if extend_window_cnf
                    [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
                end
                if(L1BS_buffer(i_surf_stacked).surface_type_flag == 4)
                    %TRP location being processed
                    disp(['Computing the results for the record' num2str(i_surf_stacked)]);
                    script_TRP_validation;
                    exit_gpp=1;
                    break;
                end
                fprintf('step(3) L1BS to L1B processing...\n')
                disp('filesid in L1B_processing'); disp(files.fid)
                fprintf('Multi-looking...\n')
                [L1BS_buffer(i_surf_stacked),L1B]   = multilooking          (L1BS_buffer(i_surf_stacked));
                if include_wfms_aligned
                    %surface alignment using window delay reference of first
                    %surface
                    fprintf('Perform surface alignment using window delay reference of first surface...\n')
                    if i_surf_stacked== 1
                        %reference surface
                        win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
                        alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;
                    end
                    [L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
                end
                fprintf('Sigma-0 scaling factor...\n')
                [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B);
                
                %plot track in the GUI
                if(files.options.plotting_flag(3))
                    plotm(L1BS_buffer(i_surf_stacked).lat_surf,L1BS_buffer(i_surf_stacked).lon_surf,'ro','Parent',options.axes);
                    h_map_points =plotm([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf],[L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lon_sar_surf],'y.','Parent',options.axes);
                    h_points = plot([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).alt_sar_surf], [L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf],'y.','Parent',options.wd_axes);
                    set(options.wd_axes, 'YLim',[min([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf])-0.05 max([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf])+0.05]);
                    pause(0.01);
                    delete(h_map_points);
                    delete(h_points);
                    plotm(L1BS_buffer(i_surf_stacked).lat_surf,L1BS_buffer(i_surf_stacked).lon_surf,'go','Parent',options.axes)
                    plot(L1BS_buffer(i_surf_stacked).alt_surf, L1BS_buffer(i_surf_stacked).lat_surf,'go','Parent',options.wd_axes);
                end
                if(options.GUI_flag)
                    drawnow();
                    stop_state = get(options.pause_button, 'Value');
                    if stop_state
                        break;
                    end
                end
                %added EM: 30.10.2016
                %% ---------------- ACDC APPLICATION ----------------------
                if ACDC_application_cnf
                    fprintf('ACDC processing...\n')
                    if i_surf_stacked==1
                        %% ----------------- FITTING PARAMS INITIALIZE ----------------------------
                        %--------------------------------------------------------------------------
                        fprintf('i_surf_stacked=1, initializing fitting parameters...\n')
                        if cnf_p_ACDC.rou_flag
                            fprintf('(MSS fitting)\n')
                            if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                                fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou cnf_p_ACDC.ini_Pu];
                            else
                                fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou];
                            end
                        else
                            fprintf('(SWH fitting)\n')
                            if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                                fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4 cnf_p_ACDC.ini_Pu];
                            else
                                fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4];
                            end
                        end
                        fit_params_ini_ACDC=[];
                        if cnf_p_ACDC.ini_Hs_rou_sliding_win_opt==1
                            accumulated_sigmaz_conv=0;
                            accumulated_SWH_ACDC=0;
                            accumulated_SSH_ACDC=0;
                        end
                        i_surf_fitted=0;
                        %% ----------------- INITIALIZE RANGE_INDEXATION --------------------------
                        %--------------------------------------------------------------------------
                        switch cnf_p_ACDC.range_index_method
                            %N_samples corresponds to non-zero padded ones
                            case 'conventional'
                                range_index=1:(N_samples*zp_fact_range_cnf);
                                delta_range=c_cst/(2.0*bw_ku_chd*zp_fact_range_cnf);
                            case 'resampling'
                                range_index=interp((1:N_samples),zp_fact_range_cnf);
                                delta_range=c_cst/(2.0*bw_ku_chd);
                        end
                        %% ----------------- LOAD THE LUTs FUNC_F0/F1 ------------------------------
                        %--------------------------------------------------------------------------
                        if (cnf_p_ACDC.lut_flag)
                            fprintf('Loading LUTs functions f0 and f1...\n')
                            %-------------- load func_f0 --------------------------------------
                            load('./inputs/LUT_f0.mat','func_f0')
                            switch cnf_p_ACDC.power_wfm_model
                                case 'complete'
                                    load('./inputs/LUT_f1.mat','func_f1')
                                case 'simple'
                                    %
                                otherwise
                                    error('No valid power waveform model')
                            end
                        end
                        geo_param_buffer_conv=[];
                        geo_param_buffer_ACDC=[];
                    end
                    fprintf('ACDC stack processing...\n')
                    [L1B,geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,i_surf_fitted,flag_exit]  = ACDC_stack_bis (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B,...
                        range_index,delta_range,i_surf_stacked,...
                        geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,...
                        accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,i_surf_fitted,...
                        func_f0,func_f1,files);
                    
                    Hs_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.Hs_conv;
                    Hs_ACDC(i_surf_stacked)=L1B.ACDC.Hs;
                    SSH_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.SSH_conv;
                    SSH_ACDC(i_surf_stacked)=L1B.ACDC.SSH;
                end
                
                %% Writting routines
                if(options.writting_flag(2))
                    if i_surf_stacked==1
                        files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
                    end
                    prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
                end
                if(options.writting_flag(1))
                    if i_surf_stacked==1
                        files.ncid_BS = netcdf.open(files.filename_netCDF_BS,'WRITE');
                    end
                    prepare_NetCDF_L1BS_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
                end
                
                indexes_to_remove = find(burst_indexes_in_memory(1:end) < min(L1BS_buffer(i_surf_stacked).burst_index));
                
                if(~isempty(indexes_to_remove))
                    %[L1A_buffer(burst_indexes_in_memory(indexes_to_remove))] = empty_L1A_struct(L1A_buffer(burst_indexes_in_memory(indexes_to_remove)));
                    %burst_indexes_in_memory(indexes_to_remove)=[];
                    
                end
                %[L1BS_buffer(i_surf_stacked)] = empty_L1BS_struct(L1BS_buffer(i_surf_stacked));
                i_surf_stacked = i_surf_stacked +1;
                if(i_surf_stacked>N_surfs_loc_estimated)
                    disp('exit process');
                    exit_gpp=1;
                    break;
                end
                
            end
        end
        
    end
    i_burst=i_burst+1;
end

last_burst=i_burst_focussed;
% Last bursts and stacks
if(~exit_gpp)
    fprintf('Processing last bursts and stacks...\n')
    for i_burst_focussed = last_burst:N_bursts
        %progressbar([],[],i_burst_focussed/N_bursts,[]);
        [L1A_buffer(i_burst_focussed)]      = beam_angles (L1A_buffer(i_burst_focussed),   L1BS_buffer, N_total_surf_loc,i_surf_stacked);
        [L1A_buffer(i_burst_focussed)]      = azimuth_processing    (L1A_buffer(i_burst_focussed));
        if(L1A_buffer(i_burst_focussed).surf_loc_index(1)~=i_surf_stacked)
            %             if(i_burst_focussed > (N_bursts_cycle_chd*N_ku_pulses_burst_chd))
            %progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
            [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
            [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
            [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
            if extend_window_cnf
                [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
            end
            if(L1BS_buffer(i_surf_stacked).surface_type_flag == 4)
                %TRP location being processed
                disp(['Computing the results for the record' num2str(i_surf_stacked)]);
                script_TRP_validation;
                exit_gpp=1;
                break;
            end
            [L1BS_buffer(i_surf_stacked),L1B]   = multilooking          (L1BS_buffer(i_surf_stacked));
            if include_wfms_aligned
                if i_surf_stacked== 1
                    %reference surface
                    win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
                    alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;
                end
                [L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
            end
            [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B);
            if(files.options.plotting_flag(3))
                plotm(L1BS_buffer(i_surf_stacked).lat_surf,L1BS_buffer(i_surf_stacked).lon_surf,'ro','Parent',options.axes);
                h_map_points =plotm([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf],[L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lon_sar_surf],'y.','Parent',options.axes);
                h_points = plot([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).alt_sar_surf], [L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf],'y.','Parent',options.wd_axes);
                set(options.wd_axes, 'YLim',[min([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf])-0.05 max([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf])+0.05]);
                pause(0.01);
                delete(h_map_points);
                delete(h_points);
                plotm(L1BS_buffer(i_surf_stacked).lat_surf,L1BS_buffer(i_surf_stacked).lon_surf,'go','Parent',options.axes)
                plot(L1BS_buffer(i_surf_stacked).alt_surf, L1BS_buffer(i_surf_stacked).lat_surf,'go','Parent',options.wd_axes);
            end
            if(options.GUI_flag)
                drawnow();
                stop_state = get(options.pause_button, 'Value');
                if stop_state
                    break;
                end
            end
            %added EM: 30.10.2016
            %% ---------------- ACDC APPLICATION ----------------------
            if ACDC_application_cnf
                if i_surf_stacked==1
                    %% ----------------- FITTING PARAMS INITIALIZE ----------------------------
                    %--------------------------------------------------------------------------
                    if cnf_p_ACDC.rou_flag
                        if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou cnf_p_ACDC.ini_Pu];
                        else
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou];
                        end
                    else
                        if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4 cnf_p_ACDC.ini_Pu];
                        else
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4];
                        end
                    end
                    fit_params_ini_ACDC=[];
                    if cnf_p_ACDC.ini_Hs_rou_sliding_win_opt==1
                        accumulated_sigmaz_conv=0;
                        accumulated_SWH_ACDC=0;
                        accumulated_SSH_ACDC=0;
                    end
                    i_surf_fitted=0;
                    %% ----------------- INITIALIZE RANGE_INDEXATION --------------------------
                    %--------------------------------------------------------------------------
                    switch cnf_p_ACDC.range_index_method
                        %N_samples corresponds to non-zero padded ones
                        case 'conventional'
                            range_index=1:(N_samples*zp_fact_range_cnf);
                            delta_range=c_cst/(2.0*bw_ku_chd*zp_fact_range_cnf);
                        case 'resampling'
                            range_index=interp((1:N_samples),zp_fact_range_cnf);
                            delta_range=c_cst/(2.0*bw_ku_chd);
                    end
                    %% ----------------- LOAD THE LUTs FUNC_F0/F1 ------------------------------
                    %--------------------------------------------------------------------------
                    if (cnf_p_ACDC.lut_flag)
                        %-------------- load func_f0 --------------------------------------
                        load('./inputs/LUT_f0.mat','func_f0')
                        switch cnf_p_ACDC.power_wfm_model
                            case 'complete'
                                load('./inputs/LUT_f1.mat','func_f1')
                            case 'simple'
                                %
                            otherwise
                                error('No valid power waveform model')
                        end
                    end
                    geo_param_buffer_conv=[];
                    geo_param_buffer_ACDC=[];
                end
                [L1B,geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,i_surf_fitted,flag_exit]  = ACDC_stack_bis (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B,...
                    range_index,delta_range,i_surf_stacked,...
                    geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,...
                    accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,i_surf_fitted,...
                    func_f0,func_f1,files);
                
                Hs_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.Hs_conv;
                Hs_ACDC(i_surf_stacked)=L1B.ACDC.Hs;
                SSH_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.SSH_conv;
                SSH_ACDC(i_surf_stacked)=L1B.ACDC.SSH;
            end
            
            %% Writting routines
            if(options.writting_flag(2))
                if i_surf_stacked==1
                    files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
                end
                prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
            end
            if(options.writting_flag(1))
                if i_surf_stacked==1
                    files.ncid_BS = netcdf.open(files.filename_netCDF_BS,'WRITE');
                end
                prepare_NetCDF_L1BS_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
            end
            
            indexes_to_remove = find(burst_indexes_in_memory(1:end) < min(L1BS_buffer(i_surf_stacked).burst_index));
            
            if(~isempty(indexes_to_remove))
                %[L1A_buffer(burst_indexes_in_memory(indexes_to_remove))] = empty_L1A_struct(L1A_buffer(burst_indexes_in_memory(indexes_to_remove)));
                %burst_indexes_in_memory(indexes_to_remove)=[];
            end
            %[L1BS_buffer(i_surf_stacked)] = empty_L1BS_struct(L1BS_buffer(i_surf_stacked));
            i_surf_stacked = i_surf_stacked +1;
            if(i_surf_stacked>N_surfs_loc_estimated)
                disp('exit process');
                break;
            end
            
        end
    end
    %stack the last surfaces
    last_stack = i_surf_stacked;
    % Last stacks
    if(~exit_gpp)
        for i_surf_stacked = last_stack:N_total_surf_loc
            
            %progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
                [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
                
                [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
                [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
                
                
                if extend_window_cnf
                    [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
                end
                if(L1BS_buffer(i_surf_stacked).surface_type_flag == 4&& exit_gpp==0)
                    %TRP location being processed
                    disp(['Computing the results for the record' num2str(i_surf_stacked)]);
                    script_TRP_validation;
                    exit_gpp=1;
                    break;
                end
                [L1BS_buffer(i_surf_stacked),L1B]   = multilooking          (L1BS_buffer(i_surf_stacked));
            if include_wfms_aligned
                if i_surf_stacked== 1
                    %reference surface
                    win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
                    alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;
                end
                [L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
            end
            [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B);
            if(files.options.plotting_flag(3))
                plotm(L1BS_buffer(i_surf_stacked).lat_surf,L1BS_buffer(i_surf_stacked).lon_surf,'ro','Parent',options.axes);
                h_map_points =plotm([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf],[L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lon_sar_surf],'y.','Parent',options.axes);
                h_points = plot([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).alt_sar_surf], [L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf],'y.','Parent',options.wd_axes);
                set(options.wd_axes, 'YLim',[min([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf])-0.05 max([L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(1:end)).lat_sar_surf])+0.05]);
                drawnow();
                delete(h_map_points);
                delete(h_points);
                plotm(L1BS_buffer(i_surf_stacked).lat_surf,L1BS_buffer(i_surf_stacked).lon_surf,'go','Parent',options.axes)
                plot(L1BS_buffer(i_surf_stacked).alt_surf, L1BS_buffer(i_surf_stacked).lat_surf,'go','Parent',options.wd_axes);
                
            end
            if(options.GUI_flag)
                drawnow();
                stop_state = get(options.pause_button, 'Value');
                if stop_state
                    break;
                end
            end
            %added EM: 30.10.2016
            %% ---------------- ACDC APPLICATION ----------------------
            if ACDC_application_cnf
                if i_surf_stacked==1
                    %% ----------------- FITTING PARAMS INITIALIZE ----------------------------
                    %--------------------------------------------------------------------------
                    if cnf_p_ACDC.rou_flag
                        if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou cnf_p_ACDC.ini_Pu];
                        else
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou];
                        end
                    else
                        if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4 cnf_p_ACDC.ini_Pu];
                        else
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4];
                        end
                    end
                    fit_params_ini_ACDC=[];
                    if cnf_p_ACDC.ini_Hs_rou_sliding_win_opt==1
                        accumulated_sigmaz_conv=0;
                        accumulated_SWH_ACDC=0;
                        accumulated_SSH_ACDC=0;
                    end
                    i_surf_fitted=0;
                    %% ----------------- INITIALIZE RANGE_INDEXATION --------------------------
                    %--------------------------------------------------------------------------
                    switch cnf_p_ACDC.range_index_method
                        %N_samples corresponds to non-zero padded ones
                        case 'conventional'
                            range_index=1:(N_samples*zp_fact_range_cnf);
                            delta_range=c_cst/(2.0*bw_ku_chd*zp_fact_range_cnf);
                        case 'resampling'
                            range_index=interp((1:N_samples),zp_fact_range_cnf);
                            delta_range=c_cst/(2.0*bw_ku_chd);
                    end
                    %% ----------------- LOAD THE LUTs FUNC_F0/F1 ------------------------------
                    %--------------------------------------------------------------------------
                    if (cnf_p_ACDC.lut_flag)
                        %-------------- load func_f0 --------------------------------------
                        load('./inputs/LUT_f0.mat','func_f0')
                        switch cnf_p_ACDC.power_wfm_model
                            case 'complete'
                                load('./inputs/LUT_f1.mat','func_f1')
                            case 'simple'
                                %
                            otherwise
                                error('No valid power waveform model')
                        end
                    end
                    geo_param_buffer_conv=[];
                    geo_param_buffer_ACDC=[];
                end
                [L1B,geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,i_surf_fitted,flag_exit]  = ACDC_stack_bis (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B,...
                    range_index,delta_range,i_surf_stacked,...
                    geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,...
                    accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,i_surf_fitted,...
                    func_f0,func_f1,files);
                
                Hs_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.Hs_conv;
                Hs_ACDC(i_surf_stacked)=L1B.ACDC.Hs;
                SSH_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.SSH_conv;
                SSH_ACDC(i_surf_stacked)=L1B.ACDC.SSH;
            end
            
            %% Writting routines
            if(options.writting_flag(2))
                if i_surf_stacked==1
                    files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
                end
                prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
            end
            if(options.writting_flag(1))
                if i_surf_stacked==1
                    files.ncid_BS = netcdf.open(files.filename_netCDF_BS,'WRITE');
                end
                prepare_NetCDF_L1BS_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
            end
            
            
                indexes_to_remove = find(burst_indexes_in_memory(1:end) < min(L1BS_buffer(i_surf_stacked).burst_index));
                
            if(~isempty(indexes_to_remove))
                [L1A_buffer(burst_indexes_in_memory(indexes_to_remove))] = empty_L1A_struct(L1A_buffer(burst_indexes_in_memory(indexes_to_remove)));
                burst_indexes_in_memory(indexes_to_remove)=[];
            end
            [L1BS_buffer(i_surf_stacked)] = empty_L1BS_struct(L1BS_buffer(i_surf_stacked));
            if(i_surf_stacked>N_surfs_loc_estimated)
                disp('exit process');
                break;
            end
            
            
        end
    end
end
%% ----------- Resize the netCDF file ------------------------------
if(options.writting_flag(2))
    netcdf.close(files.ncid); % close the already open netcdf
    files.filename_netCDF_2remove=files.filename_netCDF;
    [files] = create_NetCDF_L1B_S3(files,i_surf_stacked-1);
    resize_NetCDF_L1B_S3(files,i_surf_stacked-1);
    delete(files.filename_netCDF_2remove);
end
if(options.writting_flag(1))
    netcdf.close(files.ncid_BS); % close the already open netcdf
    files.filename_netCDF_2remove_BS=files.filename_netCDF_BS;
    [files] = create_NetCDF_L1BS_S3(files,i_surf_stacked-1);
    resize_NetCDF_L1BS_S3(files,i_surf_stacked-1);
    delete(files.filename_netCDF_2remove_BS);
end
if targz_option_active
    %remove the files DBL and HDR
    delete([files.inputPath filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name]); %DBL
    delete([files.inputPath strrep(filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name,'DBL','HDR')]); %HDR
end
clear files;
%     fclose all;

time = toc(t6);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
strcat('Processing L1A: ',filename_L1A)
disp(['Processing time for ',filename_L1A,': ',num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);



%% Retracker Processor
if(steps(4))
    [L2]   = retracker(L1B);
end
end