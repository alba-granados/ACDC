%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the AZIMUTH PROCESSING as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: The purpose of the azimuth processing is to steer the beams to the 
% different surface locations.
% 
% ----------------------------------------------------------
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/09/2013)
% 
% The purpose of the azimuth processing is to steer the beams to the 
% different surface locations.
%
% Version  record
% 1.0 2015/01/01 Imported code from S6
% 1.1 2016/02/18 Added zp_fact_azimut_cnf
% 1.2 2016/03/23 Added mode as a global
% 2.0 2016/04/01 Azimuth Processing for each burst individually. Stacking
% removed.
% 2.1 2016/04/10 For loop over the pulser removed in the prephase shift and
% hamming windowing
% 2.2 2016/04/16 hamm_win transpose removed to properly build the 64x128
% size
% 2.3 Hamming window computation removed and loaded as a global variable to 
% reduce computational load (avoid compute window each burst)
% 2.4 Changed N_samples_sar_chd N_samples
% 2.5 2016/11/16 case S3 without the fftshift?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L1A] = azimuth_processing (L1A)

    global mode mission
    global N_ku_pulses_burst_chd wv_length_ku N_samples
    global sigma_alt_surf_th_cnf force_exact_method_cnf
    global pi_cst hamming_window_cnf zp_fact_azimut_cnf hamm_win
    
    
    % Pre-allocating memory                                              
    wfm_cal_corrected_azw   = zeros (N_ku_pulses_burst_chd, N_samples);                                              
    wfm_beams_focused       = zeros (N_ku_pulses_burst_chd*zp_fact_azimut_cnf, N_samples);
    L1A.beams_focused_shifted   = zeros (N_ku_pulses_burst_chd*zp_fact_azimut_cnf, N_samples);
    
    
    if(strcmp(mode,'SIN'))
        wfm_cal_corrected_azw_2   = zeros (N_ku_pulses_burst_chd, N_samples);                                              
        wfm_beams_focused_2       = zeros (N_ku_pulses_burst_chd*zp_fact_azimut_cnf, N_samples);
        L1A.beams_focused_shifted_2   = zeros (N_ku_pulses_burst_chd*zp_fact_azimut_cnf, N_samples);    
    end
    
    k0      = 2*pi_cst / wv_length_ku;
    
    
    
        
        vel_sat = [L1A.x_vel_sat_sar, L1A.y_vel_sat_sar, L1A.z_vel_sat_sar];
        vel_sat_norm = norm(vel_sat); clear vel_sat
  
        %% 1. AZIMUTH WEIGHTING
        if(exist('azimuth_weighting_filename_chd','var'))
            [start_angle_azw, N_weights_azw, weights_azw, angles_azw] = read_weighting(azimuth_weighting_filename_chd);
        else
            weights_azw = ones(1,N_ku_pulses_burst_chd);
        end
        %
        % Missing code to interpolate the weighting to apply it properly
        %
        
        wfm_cal_corrected_azw(:,:) = L1A.wfm_cal_gain_corrected(:,:);
        if(strcmp(mode,'SIN'))
            wfm_cal_corrected_azw_2(:,:) = L1A.wfm_cal_gain_corrected_2(:,:);    
        end
        
        %% 2. METHOD SELECTION
        % Determine which method is going to be used: exact or approximate.
        % --------------------------------------------------------------------
        %  (This can be forced in the cnf file).


%         mean_alt = mean(alt_surf(L1A.surf_loc_index(1) : L1A.surf_loc_index(L1A.N_beams_sar_ku)));
%         aux=0;
%         for i_beam = 1:L1A.N_beams_sar_ku
%             aux = aux + (alt_surf(L1A.surf_loc_index(i_beam)) - mean_alt)^2;
%         end
% 
%         sigma_alt_surf = sqrt (1/L1A.N_beams_sar_ku * aux);
        sigma_alt_surf = 0;
%         force_exact_method_cnf = 1;
        if sigma_alt_surf >= sigma_alt_surf_th_cnf
            method = 1; %Exact method
        elseif(force_exact_method_cnf)
            method = 1; %Exact method
        else
            method = 0; %Approximate method
        end
        
        
        if method == 1
            
%             L1A.beam_ang_index = N_ku_pulses_burst_chd/2*zp_fact_azimut_cnf+1;
            
            %% -2.1. EXACT METHOD
            
            for i_beam = 1:L1A.N_beams_sar_ku

                switch mission
                    case 'CR2'
                
                       % Pre-azimuth FFT phase shift
                        wfm_pulses_received_shifted = L1A.wfm_cal_gain_corrected .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(i_beam))* L1A.pri_sar * (0:N_ku_pulses_burst_chd-1)).'*ones(1,N_samples)) ;

                    case {'S3_','S3A','S3B'} 
                        % Pre-azimuth FFT phase shift
                        wfm_pulses_received_shifted = L1A.wfm_cal_gain_corrected .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(i_beam))* L1A.pri_sar * (0:N_ku_pulses_burst_chd-1)).'*ones(1,N_samples)) ;

                end
                
                if(hamming_window_cnf)
                    % Hamming
                    
                    (:,:) = wfm_pulses_received_shifted(:,:) .* (hamm_win * ones(1,N_samples));

                end
                % FFT
                wfm_beams_received_shifted = fft(wfm_pulses_received_shifted,N_ku_pulses_burst_chd*zp_fact_azimut_cnf)/sqrt(N_ku_pulses_burst_chd*zp_fact_azimut_cnf);
                switch mission
                    case 'CR2'
                
                        wfm_beams_received_shifted = fftshift(wfm_beams_received_shifted,1);
                    case {'S3_','S3A','S3B'} 
                        %don't need to do any swap
                        wfm_beams_received_shifted = fftshift(wfm_beams_received_shifted,1);
                end
            
                % Select the focused beam (the central one) and store it
                wfm_beams_focused (i_beam, :) = wfm_beams_received_shifted(L1A.beam_ang_index,:);
                
                if(strcmp(mode,'SIN'))
                    
                    wfm_pulses_received_shifted_2 = L1A.wfm_cal_gain_corrected_2 .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(i_beam))* L1A.pri_sar * (0:N_ku_pulses_burst_chd-1)).'*ones(1,N_samples)) ;

                    if(hamming_window_cnf)
                    % Hamming antenna 2
                        wfm_pulses_received_shifted_2 = wfm_pulses_received_shifted_2 .* (hamm_win * ones(1,N_samples));
                    end
                    % FFT 
                    wfm_beams_received_shifted_2 = fft(wfm_pulses_received_shifted_2,N_ku_pulses_burst_chd*zp_fact_azimut_cnf)/sqrt(N_ku_pulses_burst_chd*zp_fact_azimut_cnf);
                    wfm_beams_received_shifted_2 = (wfm_beams_received_shifted_2);
                    % Select the focused beam (the central one) and store it
                    wfm_beams_focused_2 (i_beam, :) = wfm_beams_received_shifted_2(L1A.beam_ang_index,:);
                    
                end
            end
            % % Post-azimuth FFT phase shift
            % for i_beam = 1:L1A.N_beams_sar_ku
            %     L1A.beams_focused_shifted (i_beam,:) = wfm_beams_focused (i_beam,:) * exp (1i * (k0 * norm(vel_sat) * cos(L1A.beam_ang(i_surf,i_beam)) * (N_ku_pulses_burst_chd-1) * L1A.pri_sar));
            % end
            L1A.beams_focused_shifted(:,:) = wfm_beams_focused(:,:);
            if(strcmp(mode,'SIN'))
                L1A.beams_focused_shifted_2(:,:) = wfm_beams_focused_2(:,:);
            end
%             plot_wfm(:,:)=L1A.beams_focused_shifted_surf(:,:);
%             imagesc(abs(fftshift(fft(plot_wfm.'),1).'));

        elseif method == 0
            %% -2.2. APPROXIMATE METHOD

%                L1A.beam_ang_index = (N_ku_pulses_burst_chd/2)+1; 
% %               Fix azimuth processing at the end of the track 
               while(L1A.beam_ang( L1A.beam_ang_index)==0||L1A.beam_ang( L1A.beam_ang_index)==Inf)
                  L1A.beam_ang_index= L1A.beam_ang_index-1;
               end
               switch mission
                   case 'CR2'
                       
                       wfm_pulses_received_shifted = L1A.wfm_cal_gain_corrected .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(L1A.beam_ang_index))* L1A.pri_sar * (0:N_ku_pulses_burst_chd-1)).'*ones(1,N_samples)) ;
                       
                   case {'S3_','S3A','S3B'}
                       %don't need to do any swap
                       wfm_pulses_received_shifted = L1A.wfm_cal_gain_corrected .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(L1A.beam_ang_index))* L1A.pri_sar * (0:N_ku_pulses_burst_chd-1)).'*ones(1,N_samples)) ;
                       
               end
            
               
               % Pre-allocating memory
             
            
            if(hamming_window_cnf)
                % Hamming
                wfm_pulses_received_shifted(:,:) = wfm_pulses_received_shifted(:,:) .* (hamm_win * ones(1,N_samples));
                
            end
            % FFT process
            wfm_beams_focused = fft(wfm_pulses_received_shifted,N_ku_pulses_burst_chd*zp_fact_azimut_cnf)/sqrt(N_ku_pulses_burst_chd*zp_fact_azimut_cnf);
            switch mission
                case 'CR2'
                    L1A.beams_focused_shifted(:,:) = fftshift(wfm_beams_focused,1);
                        
                case {'S3_','S3A','S3B'} 
                    L1A.beams_focused_shifted(:,:) = fftshift(wfm_beams_focused,1);

                    %don't need to do any swap
             end
            
            
            if(strcmp(mode,'SIN'))
                wfm_pulses_received_shifted_2  = L1A.wfm_cal_gain_corrected_2 .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(L1A.beam_ang_index))* L1A.pri_sar * (0:N_ku_pulses_burst_chd-1)).'*ones(1,N_samples));
                
                if(hamming_window_cnf)
                    % Hamming antenna 2
                    wfm_pulses_received_shifted_2 = wfm_pulses_received_shifted_2 .* (hamm_win * ones(1,N_samples));
                end
                
                % FFT process
                wfm_beams_focused_2 = fft(wfm_pulses_received_shifted_2,N_ku_pulses_burst_chd*zp_fact_azimut_cnf)/sqrt(N_ku_pulses_burst_chd*zp_fact_azimut_cnf);
                L1A.beams_focused_shifted_2(:,:) = (fftshift(wfm_beams_focused_2,1));

            end
            
        end % if-method
    
    
%     L1A = rmfield(L1A,'wfm_cal_gain_corrected');
%     if(strcmp(mode,'SIN'))
%        	L1A = rmfield(L1A,'wfm_cal_gain_corrected_2');
%     end
   
    
    
end












