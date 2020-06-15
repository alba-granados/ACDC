
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
% v2.7 2016/11/21   writing_flags = [writing_L1B_pLRM writing_L1BS writing_L1B]

function L1B_processing(filesBulk,i_fileL1A_input,targz_option_active, writing_flags)



    global PLRM_mode
    global N_bursts_cycle_chd N_ku_pulses_burst_chd
    global extend_window_cnf include_wfms_aligned 
    global optional_ext_file_flag file_ext_string
    optional_ext_file_flag = 0;
    file_ext_string = '_a_la_CS2';
    
    exit_gpp = 0;
    PLRM_mode = 1;
    i_RC = 0;
    t6 = tic;
%% Define the files structure for the specific FBR to be processed
    % including all other files (cnf,cst,chd,....) not DBL nor HDR related to the specific
    % folder (in order no to change the functions)
    files.inputPath =filesBulk.inputPath;
    files.resultPath =filesBulk.resultPath;
    files.inputFiles =filesBulk.inputFiles(~(filesBulk.filterDATAFILES));
    if targz_option_active
        %untar the file
        untar([files.inputPath filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name],files.inputPath);
        filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name=strrep(filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name,'TGZ','DBL');
    end
    %add the current L1A to be processed
    files.inputFiles(end+1) = filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input));
    files.indexaDirs      =   find(([files.inputFiles.isdir]));
    files.indexFiles      =   find(not([files.inputFiles.isdir]));
    files.nFiles          =   length(files.indexFiles);
    
    filename_L1A=filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name;
    disp(strcat('Processing L1A: ',filename_L1A))
    disp('....')
    
    %progressbar('Burst processing','Surfaces computation','Burst focussing','Stack processing');
    i_surf              = 1;
    i_burst_focused     = 1;
    i_surf_stacked      = 1;
    % previous_gap_flag   = 0; % indicates that the burst has a gap unable to be fullfilled
 %% Read inputs   
    [steps,meteo,files,first_burst,final_burst] = read_inputs(files);
    N_bursts    = final_burst-first_burst+1;
	if(N_bursts < N_bursts_cycle_chd*N_ku_pulses_burst_chd)
		disp(['Not enough records inside mask to create a complete stack for the file ' files.filename_L1A])
		return;
	end
    original_burst_index    =   first_burst:final_burst;
    burst_indexes_in_memory = 1: N_bursts;
    N_surfs_loc_estimated = floor(N_bursts/N_bursts_cycle_chd*2.0);
    if(writing_flags(3))
		[files] = create_NetCDF_L1B_S3(files,N_surfs_loc_estimated);
    end
    if(writing_flags(3))
		[files] = create_NetCDF_L1BS_S3(files,N_surfs_loc_estimated);   
    end
    if writing_flags(3)
        [files] = create_NetCDF_L1B_pLRM(files,floor(N_bursts/N_bursts_cycle_sar_chd));
    end
    
    [L1A]           = create_L1A_struct;
    [L1A_buffer]    = create_L1A_struct;
    [L1BS_buffer]   = create_L1BS_struct;
    L1B             = [];
    
    
    
    i_burst=1;
    N_bursts_original=N_bursts;
    while(i_burst<=N_bursts)
        %progressbar(i_burst/N_bursts,[],[],[]);
        %% BURST Processor
        if (steps(1))
            [ISP]  = read_ISP_record(files,first_burst,original_burst_index(i_burst),N_bursts); %To be build
            [L1A]  = pre_datation    (ISP);
            [L1A]  = pre_win_delay   (L1A,ISP);
            [L1A]  = final_burst_datation   (L1A,ISP);
            [L1A]  = onboard_reversion   (L1A,ISP);
            [L1A]  = instrument_gain_corr   (ISP,L1A);
            [L1A]  = waveforms_correction   (ISP,L1A);
        else
            
            [L1A,files]  = read_L1A_record(files,L1A,original_burst_index(i_burst),i_burst);
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
        
        if(steps(2))
            [L1BS_buffer, out_surf] = surface_locations (L1A_buffer,L1BS_buffer, i_burst, i_surf);

            
            if(i_surf < out_surf)
                %progressbar([],i_surf/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd,[],[]);
                i_surf = out_surf;
                N_total_surf_loc = i_surf-1;
                new_surface=1;
                %             disp([i_burst N_total_surf_loc]);
            end
            
            if(i_surf > 64)
                
                if(i_burst_focused==1)% handle final bursts
                    burst_margin = i_burst;
                end
                %progressbar([],[],i_burst_focused/N_bursts,[]);
                
                [L1A_buffer(i_burst_focused)]      = beam_angles (L1A_buffer(i_burst_focused),   L1BS_buffer, N_total_surf_loc,i_surf_stacked);
                [L1A_buffer(i_burst_focused)]      = azimuth_processing    (L1A_buffer(i_burst_focused));
                i_burst_focused = i_burst_focused+1;
                
                if(L1A_buffer(i_burst_focussed-1).surf_loc_index(1)~=i_surf_stacked)
                    %             if(i_burst_focussed > (N_bursts_cycle_chd*N_ku_pulses_burst_chd))
                    %progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
                    [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
                    [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
                    [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
                    if extend_window_cnf
                        [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
                    end
                    [L1BS_buffer(i_surf_stacked),L1B]   = multilooking          (L1BS_buffer(i_surf_stacked));
                    if include_wfms_aligned
                        %surface alignment using window delay reference of first
                        %surface
                        if i_surf_stacked== 1
                            %reference surface
                            win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
                            alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;
                        end
                        [L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
                    end
                    [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1BS_buffer(i_surf_stacked),L1B);
                                       
                    
                    %% Writting routines
					
					if(writing_flags(3))
                        if i_surf_stacked==1
                            files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
                        end
                        prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
					end
					if(writing_flags(3))
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
                    i_surf_stacked = i_surf_stacked +1;
                    if(i_surf_stacked>N_surfs_loc_estimated)
                        disp('exit process');
                        exit_gpp=1;
                        break;
                    end
                    
                end
            end
            
          
        elseif L1A_buffer(i_burst).burst_sar_ku_fbr == 4
            %% pLRM chain  
            i_RC = i_RC + 1;
            L1B = plrm_chain(L1A_buffer(i_burst-3:i_burst));
            if writing_flags(3)
                if i_RC == 1
                    files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
                end
                prepare_NetCDF_L1B_pLRM (files,L1B,i_RC);
            end
            if i_burst == N_bursts
                disp('exit process');
                exit_gpp = 1;
                break;
            end
        end
    i_burst=i_burst+1;
    end
        
    %%
    %close the input binary file
    fclose(files.fid);
    last_burst=i_burst_focused;
    % Last bursts and stacks
    if(~exit_gpp)
        for i_burst_focused = last_burst:N_bursts            
            %progressbar([],[],i_burst_focused/N_bursts,[]);
            [L1A_buffer(i_burst_focused)]      = beam_angles (L1A_buffer(i_burst_focused),   L1BS_buffer, N_total_surf_loc,i_surf_stacked);
            [L1A_buffer(i_burst_focused)]      = azimuth_processing    (L1A_buffer(i_burst_focused));
            if(L1A_buffer(i_burst_focussed).surf_loc_index(1)~=i_surf_stacked)
                %             if(i_burst_focussed > (N_bursts_cycle_chd*N_ku_pulses_burst_chd))
                %progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
                [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
                [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
                [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
                if extend_window_cnf
                    [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
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
                [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1BS_buffer(i_surf_stacked),L1B);
                %% Writting routines
					
				if(writing_flags(3))
					if i_surf_stacked==1
						files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
					end
					prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
				end
				if(writing_flags(3))
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
        for i_surf_stacked = last_stack:N_total_surf_loc
                        
            %progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
            [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
            [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
            [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
            if extend_window_cnf
                [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
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
            [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1BS_buffer(i_surf_stacked),L1B);
            %% Writting routines
			if(writing_flags(3))
				if i_surf_stacked==1
					files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
				end
				prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
			end
			if(writing_flags(3))
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
%             i_surf_stacked = i_surf_stacked +1;
            if(i_surf_stacked>N_surfs_loc_estimated)
                disp('exit process');
                break;
            end
            
            
        end
        
    end
    %% ----------- Resize the netCDF file ------------------------------
    if(writing_flags(3))
        netcdf.close(files.ncid); % close the already open netcdf
        files.filename_netCDF_2remove=files.filename_netCDF;
        [files] = create_NetCDF_L1B_S3(files,i_surf_stacked-1);
        resize_NetCDF_L1B_S3(files,i_surf_stacked-1);
        delete(files.filename_netCDF_2remove);
    end
    if(writing_flags(3))
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
    fclose all;
    close all;
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