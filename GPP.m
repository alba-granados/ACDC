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
% Authors:  Albert Garcia-Mondejar / isardSAT
%           Roger Escola Jane / isardSAT
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clear global
set_default_plot
writting = 1;
exit_gpp = 0;
global zp_fact_azimut_cnf N_bursts_cycle_chd trp_flag_cnf

t6 = tic;

% define Paths
files.inputPath       =   '..\inputs\FBR_ESA\ROI_4\05\';%'..\inputs\FBR_ESA\ROI_1\10_subset\';
files.resultPath      =   '..\results\L1B_ISD\ROI_4\05\';

files.inputFiles      =   dir(files.inputPath);
files.indexaDirs      =   find(([files.inputFiles.isdir]));
files.indexFiles      =   find(not([files.inputFiles.isdir]));
files.nFiles          =   length(files.indexFiles);             % number of input files


progressbar('Burst processing','Surfaces computation','Burst focussing','Stack processing');
i_surf              = 1;
i_burst_focussed    = 1;
i_surf_stacked      = 1;
% previous_gap_flag   = 0; % indicates that the burst has a gap unable to be fullfilled 

[steps,meteo,files,first_burst,final_burst] = read_inputs(files);
N_bursts    = final_burst-first_burst+1;
original_burst_index    =   first_burst:final_burst;
burst_indexes_in_memory = 1: N_bursts;
N_surfs_loc_estimated = floor(N_bursts/N_bursts_cycle_chd)+100;
if(writting)
    [files] = create_NetCDF_L1B_S3(files,N_surfs_loc_estimated);
end
[L1A]           = create_L1A_struct;
[L1A_buffer]    = create_L1A_struct;
[L1BS_buffer]   = create_L1BS_struct;
L1B             = [];



for i_burst = 1: N_bursts
    progressbar(i_burst/N_bursts,[],[],[]);
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
        [L1A]  = read_L1A_record(files,L1A,original_burst_index(i_burst),N_bursts);
        if(i_burst > 1)
            %if there is are missing burst records in the FBR file, jumps
            %in time, they are created to provide continuity to the stacks.
            [L1A,original_burst_index,N_bursts] = check_continuity (L1A_buffer(end),L1A,original_burst_index,i_burst,N_bursts);
        end
        
        
        
        if(L1A.confi_block_degraded==1&&L1A.days==0) % handle degraded bursts
           [L1A] = fill_empty_burst(files, L1A_buffer(end),L1A,original_burst_index(i_burst),N_bursts); 
        end
    end
    L1A_buffer(i_burst) = L1A;

    if(steps(2))
        [L1BS_buffer, out_surf] = surface_locations (L1A_buffer,L1BS_buffer, i_burst, i_surf);

        if(i_surf < out_surf)
            progressbar([],i_surf/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd,[],[]);
            i_surf = out_surf;
            N_total_surf_loc = i_surf-1;
            new_surface=1;
%             disp([i_burst N_total_surf_loc]);
        end

        if(i_surf > 64)

            if(i_burst_focussed==1)% handle final bursts
                burst_margin = i_burst;
            end
            progressbar([],[],i_burst_focussed/N_bursts,[]);

            [L1A_buffer(i_burst_focussed)]      = beam_angles (L1A_buffer(i_burst_focussed),   L1BS_buffer, N_total_surf_loc,i_surf_stacked);
            [L1A_buffer(i_burst_focussed)]      = azimuth_processing    (L1A_buffer(i_burst_focussed));
            i_burst_focussed = i_burst_focussed+1;

            if(L1A_buffer(i_burst_focussed-1).surf_loc_index(1)~=i_surf_stacked)
%             if(i_burst_focussed > (N_bursts_cycle_chd*N_ku_pulses_burst_chd))
                progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
                [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
				if(L1BS_buffer(i_surf_stacked).surface_type_flag==4 || trp_flag_cnf == 0)
            
					[L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
					[L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
					[L1BS_buffer(i_surf_stacked),L1B]   = multilooking          (L1BS_buffer(i_surf_stacked)); 
					%surface alignment using window delay reference of first
					%surface
					if i_surf_stacked== 1 
						%reference surface
						win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
						alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;                                     
					end
					[L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
					[L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1BS_buffer(i_surf_stacked),L1B);
					%% Writting routines
					if(writting)
	%                     prepare_NetCDF_L1BS_S3(L1BS_buffer(i_surf_stacked));
						prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
					end
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
       
    end

end
last_burst=i_burst_focussed;
% Last bursts and stacks
if(~exit_gpp)
    for i_burst_focussed = last_burst:N_bursts
        progressbar([],[],i_burst_focussed/N_bursts,[]);
        [L1A_buffer(i_burst_focussed)]      = beam_angles (L1A_buffer(i_burst_focussed),   L1BS_buffer, N_total_surf_loc,i_surf_stacked);
        [L1A_buffer(i_burst_focussed)]      = azimuth_processing    (L1A_buffer(i_burst_focussed));

        if(L1A_buffer(i_burst_focussed).surf_loc_index(1)~=i_surf_stacked)
            progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
            [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
            if(L1BS_buffer(i_surf_stacked).surface_type_flag ==4 || trp_flag_cnf == 0)
				[L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
				[L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
				[L1BS_buffer(i_surf_stacked),L1B]   = multilooking          (L1BS_buffer(i_surf_stacked)); 
				[L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
				[L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1BS_buffer(i_surf_stacked),L1B);
				%% Writting routines
				if(writting)
		%             prepare_NetCDF_L1BS_S3(L1BS_buffer(i_surf_stacked));
					prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
				end
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
end
%% ----------- Resize the netCDF file ------------------------------
if(writting)
    files.filename_netCDF_2remove=files.filename_netCDF;
    [files] = create_NetCDF_L1B_S3(files,i_surf_stacked-1);
    resize_NetCDF_L1B_S3(files,i_surf_stacked-1);
    delete(files.filename_netCDF_2remove);
end

time = toc(t6);
minutes_reading = floor(time/60);
secs_reading = time - minutes_reading*60;
disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds processing L1B']);



%% Retracker Processor
if(steps(4))
    [L2]   = retracker(L1B);
end    
    
    
    
    


