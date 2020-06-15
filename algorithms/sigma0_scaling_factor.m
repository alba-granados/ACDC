%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 9.11 WAVEFORMS SCALING FACTOR %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of the waveforms scaling factor is to provide the L2
% processing with the waveform scaling factor in order to compute the
% backscatter coefficient of the surface from which the echo is reflected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v1.1 added S3 case. Sigma0 computed from the one comming from L1A. See
% "adaptnetCDF2internal.m"


%____________________________For each stack________________________________

function [L1BS , L1B]   = sigma0_scaling_factor (L1A, L1BS, L1B)


global N_ku_pulses_burst_chd
global power_tx_ant_ku_chd antenna_gain_ku_chd wv_length_ku
global pulse_length_chd chirp_slope_ku_chd
global c_cst pi_cst
global hamming_window_cnf
global earth_radius_cst N_max_beams_stack_chd
global PTR_width_chd 
global mission


switch mission
    case 'CR2' 

L1BS.wfm_scaling_factor_sar_ku_beam = zeros (1,N_max_beams_stack_chd);
azimuth_distance = zeros (1,N_max_beams_stack_chd);
range_distance = zeros (1,N_max_beams_stack_chd);
surface_area = zeros (1,N_max_beams_stack_chd);
L1B.wfm_scaling_factor_sar_ku = zeros (1,1);
    
    % 1. Preliminary factors computation
    %*********************************************************
    %**** Moved to Instrument Gain Corrections algorithm *****
    %*********************************************************
%     %   > Instrument attenuation averaging
%     aux = 0;
%     for 1:L1BS.N_beams_stack = 1:L1BS.N_beams_stack
%         aux = aux + att_corr_instr_sar(L1BS.burst_index1:L1BS.N_beams_stack));
%     end
%     att_corr_instr_sar_mean = aux / L1BS.N_beams_stack;
    %*********************************************************
    %*********************************************************
    %*********************************************************
    
    %   > Surface area (of every beam)
    % Pre-allocating memory
        norm_vel_sat(1:L1BS.N_beams_stack) = sqrt(L1BS.x_vel_sat_beam(1:L1BS.N_beams_stack).^2+ L1BS.y_vel_sat_beam(1:L1BS.N_beams_stack).^2+ L1BS.z_vel_sat_beam(1:L1BS.N_beams_stack).^2);
        azimuth_distance (1:L1BS.N_beams_stack) = (1+L1BS.range_sat_surf(1:L1BS.N_beams_stack)./earth_radius_cst)*wv_length_ku .* L1BS.range_sat_surf(1:L1BS.N_beams_stack) ./ L1BS.pri_sar_sat_beam(1:L1BS.N_beams_stack) ./ (2 * (norm_vel_sat) .* N_ku_pulses_burst_chd);
%         range_distance (1:L1BS.N_beams_stack) = 2 * sqrt (c_cst * L1BS.range_sat_surf(1:L1BS.N_beams_stack) * (1/(pulse_length_chd * chirp_slope_ku_chd)).*...
%                                     earth_radius_cst./(earth_radius_cst + L1BS.range_sat_surf(1:L1BS.N_beams_stack)) );
        range_distance (1:L1BS.N_beams_stack) = 2 * sqrt (c_cst * L1BS.range_sat_surf(1:L1BS.N_beams_stack) * (PTR_width_chd).*...
                                    earth_radius_cst./(earth_radius_cst + L1BS.range_sat_surf(1:L1BS.N_beams_stack)) );
        if hamming_window_cnf
            wf=1.486*0.92;%widening factor as denoted by Dinardo
        else
            wf=1.0;
        end
        %area_0.886*wf
        surface_area(1:L1BS.N_beams_stack) = (wf.*azimuth_distance(1:L1BS.N_beams_stack) .* range_distance(1:L1BS.N_beams_stack)).*0.886;
    
    surface_area_LRM= pi_cst*(0.5* min(range_distance(find(range_distance(:)))))^2;
    
    % 3. Compute the waveform scaling factor with external contributions for each beam
    aux = 0;
%         L1BS.wfm_scaling_factor_sar_ku_beam (1:L1BS.N_beams_stack) = 10*log10(64) + 30*log10(pi_cst)...
%             + 40*log10(L1BS.range_sat_surf(1:L1BS.N_beams_stack)) - 10*log10(power_tx_ant_ku_chd) - 2*(antenna_gain_ku_chd)...
%             - 20*log10(wv_length_ku) - 10*log10(mean(surface_area(1:L1BS.N_beams_stack)))- 10*log10(64);
        
        %added by EM 30.11.2015: consider the surface and not the mean over
        %the different beams & compensate for norm. in fft range
        %compression & TBP or pulse compression gain
        L1BS.wfm_scaling_factor_sar_ku_beam (1:L1BS.N_beams_stack) = 10*log10(64) + 30*log10(pi_cst)...
            + 40*log10(L1BS.range_sat_surf(1:L1BS.N_beams_stack)) - 10*log10(power_tx_ant_ku_chd) - 2*(antenna_gain_ku_chd)...
            - 20*log10(wv_length_ku) - 10*log10(surface_area(1:L1BS.N_beams_stack));
    
    
%         L1B.wfm_scaling_factor_sar_ku = 10*log10(64) + 30*log10(pi_cst)...
%             + 40*log10(mean(L1BS.range_sat_surf(:))) - 10*log10(power_tx_ant_ku_chd) - 2*(antenna_gain_ku_chd)...
%             - 20*log10(wv_length_ku) - 10*log10(min(surface_area(:)));
%         
%         aux = aux + L1BS.wfm_scaling_factor_sar_ku_beam(1:L1BS.N_beams_stack);
    
    
    
    % 4. Compute the averaged waveform scaling factor with external contributions
    L1B.wfm_scaling_factor_sar_ku = mean(L1BS.wfm_scaling_factor_sar_ku_beam (1:L1BS.N_beams_stack));


    % 2. Compute scaling factor only with instrumental factors
%     L1B.wfm_scaling_factor_sar_ku_instr = att_corr_instr_sar_mean + 10*log10(4) + 10*log10(pi_cst)...
%         -(power_tx_ant_chd) - 2*(antenna_gain_ku_chd) - 20*log10(wv_length_ku) - (tx_rx_gain_ground_chd);
    
    L1B.wfm_scaling_factor_sar_ku_instr = 10*log10(4) + 10*log10(pi_cst)...
        -10*log10(power_tx_ant_ku_chd) - 2*(antenna_gain_ku_chd) - 20*log10(wv_length_ku)-10*log10(64) ;
    
     case 'S3_'
         
         
             
        
         L1B.wfm_scaling_factor_sar_ku = mean([L1A([L1BS.burst_index]).instrument_sigma0_correction_rx]);
         
         
         L1B.wfm_scaling_factor_sar_ku_instr = mean([L1A([L1BS.burst_index]).instrument_sigma0_correction_tx_rx]);
         
         
end
    
    
    
    
    
end


