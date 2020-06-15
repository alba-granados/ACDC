%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the STACKING PROCESSING as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% Gather all the parameters from the burst locations to the corresponding
% suface positions
%
% Version  record
% 1.0 2015/01/01 Imported code from S6
% 2.0 2016/04/01 Cutted from azimuth processing algorithm. 
% 2.1 2016/04/25 Added L1BS.surface_type_flag when no trp_flag. Surface location algorithm change this if trp_flag==1;
% 2.0 2016/05/03 Perform the burst and beam indexation of the stacking 
% using matrix notation instead of loops
% 2.1 Added L1BS.lat_sat_beam L1BS.lon_sat_beam L1BS.alt_sat_beam



function [L1BS] = stacking(L1A_buffer,L1BS)

    global pi_cst mode
    global N_ku_pulses_burst_chd zp_fact_azimut_cnf

    
    %% 1. PERFORM THE ASSOCIATION


%     B = 0;
%     for i_burst = 1:length(L1A_buffer)
%         pre_B = B;
%         %           for i_beam = L1A_buffer.start_beam(i_burst):L1A_buffer.end_beam(i_burst)
%         aux=L1A_buffer(i_burst).N_beams_sar_ku;
%         for i_beam =1:aux
%             if L1A_buffer(i_burst).surf_loc_index(i_beam) == L1BS.surf_counter
%                 B = B + 1;
%                 L1BS.burst_index(B) = i_burst;
%                 L1BS.beam_index(B) = i_beam;
%                 break;
%             end
%         end
%         if pre_B == B && B > 0 % Check if the current surface has not been illuminated yet or no longer illuminated.
%             break;
%         end
%     end
% 
%      % Determine the real number of beams per stack
%    L1BS.N_beams_stack = B;
% %     tmp_beams_surf = L1A_buffer(L1BS.burst_index(1:L1BS.N_beams_stack)).beams_focused_shifted;
% %     L1BS.beams_surf(1:L1BS.N_beams_stack,:)        = tmp_beams_surf(L1BS.beam_index(1:L1BS.N_beams_stack),:);
   
%    indx_non_zeroBeams=find([L1A_buffer(:).N_beams_sar_ku]>0);
%    surf_loc_index_matrix=reshape([L1A_buffer(indx_non_zeroBeams).surf_loc_index],[N_ku_pulses_burst_chd.*zp_fact_azimut_cnf,length(indx_non_zeroBeams)]).';
%    [L1BS.beam_index,L1BS.burst_index]= find(surf_loc_index_matrix.'==L1BS.surf_counter);
%    L1BS.burst_index=L1BS.burst_index-1+indx_non_zeroBeams(1); %indx_non_zeroBeams(1) to take into account emptied bursts
%    L1BS.N_beams_stack = length(L1BS.burst_index);
   
   indx_non_zeroBeams=[L1A_buffer(:).N_beams_sar_ku]>0;
   surf_loc_index_matrix=reshape([L1A_buffer(indx_non_zeroBeams).surf_loc_index],[N_ku_pulses_burst_chd.*zp_fact_azimut_cnf,nnz(indx_non_zeroBeams)]).';
   [L1BS.beam_index,L1BS.burst_index]= find(surf_loc_index_matrix.'==L1BS.surf_counter);
   L1BS.burst_index=L1BS.burst_index-1+find(indx_non_zeroBeams,1,'first'); %indx_non_zeroBeams(1) to take into account emptied bursts
   L1BS.N_beams_stack = length(L1BS.burst_index);
   
%     aux=reshape(permute(reshape([L1A_buffer(indx_non_zeroBeams).beams_focused_shifted],[N_ku_pulses_burst_chd.*zp_fact_azimut_cnf,N_samples_sar_chd,length(indx_non_zeroBeams)]),[3 1 2]),[N_ku_pulses_burst_chd.*zp_fact_azimut_cnf.*length(indx_non_zeroBeams),N_samples_sar_chd]);
%     L1BS.beams_surf=aux((L1BS.beam_index-1).*length(indx_non_zeroBeams)+(L1BS.burst_index-indx_non_zeroBeams(1)+1),:);
%     % matrix way to create the stack 
%     aux=reshape([L1A_buffer(indx_non_zeroBeams).beams_focused_shifted],[N_ku_pulses_burst_chd.*zp_fact_azimut_cnf,N_samples_sar_chd,length(indx_non_zeroBeams)]);  
%     indexes=sub2ind(size(aux),...
%         reshape(ones(N_samples_sar_chd,1)*L1BS.beam_index.',[1,N_samples_sar_chd*L1BS.N_beams_stack]),...
%         reshape((1:N_samples_sar_chd).'*ones(1,L1BS.N_beams_stack),[1,N_samples_sar_chd*L1BS.N_beams_stack]),...
%         reshape(ones(N_samples_sar_chd,1)*(L1BS.burst_index-indx_non_zeroBeams(1)+1).',[1,N_samples_sar_chd*L1BS.N_beams_stack]));
%     L1BS.beams_surf=reshape(aux(indexes),[N_samples_sar_chd,L1BS.N_beams_stack]).';
    
    L1BS.T0_sar_surf                = [L1A_buffer(L1BS.burst_index).T0_sar];
    
    L1BS.doppler_ang_surf           = [L1A_buffer(L1BS.burst_index).doppler_ang_sar_sat];
    L1BS.ProcessID_beam             = [L1A_buffer(L1BS.burst_index).ProcessID];
    L1BS.lat_sat_beam               = [L1A_buffer(L1BS.burst_index).lat_sar_sat];
    L1BS.lon_sat_beam               = [L1A_buffer(L1BS.burst_index).lon_sar_sat];
    L1BS.alt_sat_beam               = [L1A_buffer(L1BS.burst_index).alt_sar_sat];

    L1BS.x_vel_sat_beam             = [L1A_buffer(L1BS.burst_index).x_vel_sat_sar];
    L1BS.y_vel_sat_beam             = [L1A_buffer(L1BS.burst_index).y_vel_sat_sar];
    L1BS.z_vel_sat_beam             = [L1A_buffer(L1BS.burst_index).z_vel_sat_sar];
    L1BS.x_sar_sat_beam             = [L1A_buffer(L1BS.burst_index).x_sar_sat];
    L1BS.y_sar_sat_beam             = [L1A_buffer(L1BS.burst_index).y_sar_sat];
    L1BS.z_sar_sat_beam             = [L1A_buffer(L1BS.burst_index).z_sar_sat];
    L1BS.win_delay_sar_ku_beam      = [L1A_buffer(L1BS.burst_index).win_delay_sar_ku];
    L1BS.pri_sar_sat_beam           = [L1A_buffer(L1BS.burst_index).pri_sar];
    for i_beam = 1:L1BS.N_beams_stack
        L1BS.beams_surf(i_beam,:)        = L1A_buffer(L1BS.burst_index(i_beam)).beams_focused_shifted(L1BS.beam_index(i_beam),:);
        L1BS.beam_ang_surf(i_beam)       = L1A_buffer(L1BS.burst_index(i_beam)).beam_ang(L1BS.beam_index(i_beam));
%         L1BS.T0_sar_surf(i_beam)         = L1A_buffer(L1BS.burst_index(i_beam)).T0_sar;
        L1BS.pointing_ang_surf(i_beam)   = pi_cst/2 + L1A_buffer(L1BS.burst_index(i_beam)).doppler_ang_sar_sat-L1BS.beam_ang_surf(i_beam)-L1A_buffer(L1BS.burst_index(i_beam)).pitch_sar; 
        L1BS.look_ang_surf(i_beam)       = L1BS.pointing_ang_surf(i_beam)+L1A_buffer(L1BS.burst_index(i_beam)).pitch_sar;
%         L1BS.doppler_ang_surf(i_beam)    = L1A_buffer(L1BS.burst_index(i_beam)).doppler_ang_sar_sat;
        L1BS.Gap_flag(i_beam)            = L1A_buffer(L1BS.burst_index(i_beam)).ProcessID>0;
%         L1BS.ProcessID_beam(i_beam)      = L1A_buffer(L1BS.burst_index(i_beam)).ProcessID;
        
%         L1BS.x_vel_sat_beam(i_beam)             = L1A_buffer(L1BS.burst_index(i_beam)).x_vel_sat_sar;
%         L1BS.y_vel_sat_beam(i_beam)             = L1A_buffer(L1BS.burst_index(i_beam)).y_vel_sat_sar; 
%         L1BS.z_vel_sat_beam(i_beam)             = L1A_buffer(L1BS.burst_index(i_beam)).z_vel_sat_sar;
%         L1BS.x_sar_sat_beam(i_beam)             = L1A_buffer(L1BS.burst_index(i_beam)).x_sar_sat;
%         L1BS.y_sar_sat_beam(i_beam)             = L1A_buffer(L1BS.burst_index(i_beam)).y_sar_sat; 
%         L1BS.z_sar_sat_beam(i_beam)             = L1A_buffer(L1BS.burst_index(i_beam)).z_sar_sat;
%         L1BS.win_delay_sar_ku_beam(i_beam)      = L1A_buffer(L1BS.burst_index(i_beam)).win_delay_sar_ku; 
%         L1BS.pri_sar_sat_beam(i_beam)           = L1A_buffer(L1BS.burst_index(i_beam)).pri_sar;
        
        
       
    end

    if(strcmp(mode,'SIN'))
        for i_beam = 1:L1BS.N_beams_stack
            L1BS.beams_surf_2(i_beam,:) 	= L1A_buffer(L1BS.burst_index(i_beam)).beams_focused_shifted_2(L1BS.beam_index(i_beam),:);
        end

    end
    if(isempty(L1BS.surface_type_flag))
        [~,beam_index_nadir]=min(abs(-pi/2+L1BS.beam_ang_surf'));
        burst_index_nadir= L1BS.burst_index(beam_index_nadir);
        L1BS.surface_type_flag=L1A_buffer(burst_index_nadir).surface_type_flag_bursts;
    end
    L1BS.ProcessID_surf = max([L1A_buffer(L1BS.burst_index(1:L1BS.N_beams_stack)).ProcessID]);
    %L1BS.look_ang_surf  =   L1BS.pointing_ang_surf;


    
    
end