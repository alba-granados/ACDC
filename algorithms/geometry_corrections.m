%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop
% This code implements the GEOMETRY CORRECTIONS as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: The purpose of the range compression is to perform range 
% compression of the input bursts (with an FFT) and then generate the 
% power waveforms.
% 
% ----------------------------------------------------------
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (02/06/2015)
% 
% 
% The purpose of the geometry corrections is to compute and apply:
% > Doppler correction
% > Slant range correction.
% > Window delay Misalignments
%
% Version  record
% 1.0 2015/01/01    Imported code from S6
% 1.1 2015/06/01    Added computation of N_windows.
%                   Application of shif splitted in coarse and fine to avoid wrapps   
% 1.2 2016/02/18    Added zp_fact_azimut_cnf
% 1.3 2016/03/23    Added mode as global
% 2.0 2016/04/04    Computation for each stack record
% 2.1 2016/04/10    Loops Removed
% 2.2 2016/04/25 	fixed win_delay_ref_cnf == 1 case
% 2.3 2016/05/18    Changed N_samples_sar_chd by N_samples
% 2.4 2016/06/04    Fixed Slant Range changed win_delay_ref for L1BS.win_delay_surf
% 2.5 2016/07/13    Flags apply corrections apply_doppler apply_slant apply_wd
% 2.6 2017/03/23    Added cnf flag processing_mode_cnf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [L1BS]  = geometry_corrections (L1BS)
                                                                                             
global mode processing_mode_cnf
global pulse_length_chd N_samples bw_ku_chd
global wv_length_ku zp_fact_range_cnf
global c_cst pi_cst prf_sar_nom N_total_bursts_sar_ku 
global avoid_wrapps_cnf T0_chd 
global win_delay_ref_cnf elevation_ref_cnf N_max_beams_stack_chd
global apply_doppler apply_slant apply_wd


% Pre-allocating memory
doppler_range           = zeros(1,N_max_beams_stack_chd);
i_samples               = (0:(N_samples-1));
win_delay_ref           = 0;
win_delay_ref_index     = 0;
win_delay_ref_index_2   = 0;
beam_power              = zeros(1,N_max_beams_stack_chd);
beam_power2             = zeros(1,N_max_beams_stack_chd);


% L1BS.doppler_corr       = zeros(1,max(L1BS.N_beams_stack));
% L1BS.range_sat_surf          = zeros(1,max(L1BS.N_beams_stack));
% L1BS.slant_range_corr_time   = zeros(1,max(L1BS.N_beams_stack));
% L1BS.slant_range_corr        = zeros(1,max(L1BS.N_beams_stack));
% L1BS.wd_corr                 = zeros(1,max(L1BS.N_beams_stack));
% L1BS.shift                   = zeros(1,max(L1BS.N_beams_stack));
% L1BS.shift_coarse            = zeros(1,max(L1BS.N_beams_stack));
%L1BS.good_samples            = zeros(max(L1BS.N_beams_stack),N_samples*zp_fact_range_cnf);
%              
    
   %0 = normal alignment
    %1 = optimal alignment (wrt max beam accumulated power)
    %2 = force the alignment to the closest of a given window_delay
    %3 = force to the first one.
    %4 = force to the minimum one.
   
    %% 1. BEAM POWER
    if win_delay_ref_cnf == 0
        win_delay_ref = L1BS.win_delay_surf;
		[~,L1BS.beam_ref] = min(abs(win_delay_ref-L1BS.win_delay_sar_ku_beam(1:L1BS.N_beams_stack)));
    elseif win_delay_ref_cnf ==1
        beams_surf_aux1_fft = fftshift(fft(L1BS.beams_surf(:,:).'),1).';
        
        beam_power(1:L1BS.N_beams_stack) = sum(abs(beams_surf_aux1_fft(1:L1BS.N_beams_stack,:)).^2,2);
        
        [~,L1BS.beam_ref] = max(beam_power(1:L1BS.N_beams_stack));
        win_delay_ref = L1BS.win_delay_sar_ku_beam(L1BS.beam_ref);
        L1BS.win_delay_surf = win_delay_ref;
    elseif win_delay_ref_cnf==2
       %need testing, not working
%         win_delay_ref_index_beams(:) = (L1A.alt_sar_sat((1:L1BS.N_beams_stack))- elevation_ref_cnf) * 2 / c_cst;
        [~,L1BS.beam_ref] = min(abs(win_delay_ref-L1BS.win_delay_sar_ku_beam(1:L1BS.N_beams_stack)));
    elseif win_delay_ref_cnf==3 
		L1BS.beam_ref=1;
        win_delay_ref = L1BS.win_delay_sar_ku_beam(L1BS.beam_ref);
		L1BS.win_delay_surf = win_delay_ref;
    elseif win_delay_ref_cnf==4
        [win_delay_ref,L1BS.beam_ref] = min(L1BS.win_delay_sar_ku_beam(1:L1BS.N_beams_stack));
		L1BS.win_delay_surf = win_delay_ref;
	elseif win_delay_ref_cnf==5
        [win_delay_ref,L1BS.beam_ref] = max(L1BS.win_delay_sar_ku_beam(1:L1BS.N_beams_stack));
		L1BS.win_delay_surf = win_delay_ref;
    end
      
    
        %% 2. DOPPLER CORRECTION
       norm_vel_sat(1:L1BS.N_beams_stack) = sqrt(L1BS.x_vel_sat_beam(1:L1BS.N_beams_stack).^2+ L1BS.y_vel_sat_beam(1:L1BS.N_beams_stack).^2+ L1BS.z_vel_sat_beam(1:L1BS.N_beams_stack).^2);

        
        % Range correction computation
        doppler_range(1:L1BS.N_beams_stack) = (- c_cst / wv_length_ku .* norm_vel_sat(1:L1BS.N_beams_stack) .* cos(L1BS.beam_ang_surf(1:L1BS.N_beams_stack))) .* pulse_length_chd / bw_ku_chd;
        % Doppler correction computation
        L1BS.doppler_corr(1:L1BS.N_beams_stack) = 2 / c_cst .* doppler_range(1:L1BS.N_beams_stack) ./ L1BS.T0_sar_surf(1:L1BS.N_beams_stack);
%         L1BS.doppler_corr(i_beam) = 0; %testing
        
        %% 3. SLANT RANGE CORRECTION
        % Range computation
        range_sat_surf(:,1:L1BS.N_beams_stack) = [L1BS.x_sar_sat_beam((1:L1BS.N_beams_stack)); ...
                                               L1BS.y_sar_sat_beam((1:L1BS.N_beams_stack)); ...
                                               L1BS.z_sar_sat_beam((1:L1BS.N_beams_stack))] - [L1BS.x_surf; L1BS.y_surf; L1BS.z_surf]*ones(1,L1BS.N_beams_stack);
        L1BS.range_sat_surf(1:L1BS.N_beams_stack) = sqrt((range_sat_surf(1,:)).^2+ (range_sat_surf(2,:)).^2+ (range_sat_surf(3,:)).^2);
        % Range delay computation
        L1BS.slant_range_corr_time(1:L1BS.N_beams_stack) = L1BS.win_delay_surf - (L1BS.range_sat_surf(1:L1BS.N_beams_stack) .* 2 / c_cst);
        L1BS.slant_range_corr(1:L1BS.N_beams_stack) = L1BS.slant_range_corr_time(1:L1BS.N_beams_stack) ./ L1BS.T0_sar_surf(1:L1BS.N_beams_stack);
    
        %% 4. WINDOW DELAY MISALIGNMENTS CORRECTION
        L1BS.wd_corr(1:L1BS.N_beams_stack) = - (win_delay_ref - L1BS.win_delay_sar_ku_beam((1:L1BS.N_beams_stack))) ./ L1BS.T0_sar_surf(1:L1BS.N_beams_stack);
    
    
        L1BS.shift = L1BS.doppler_corr.*apply_doppler + L1BS.slant_range_corr.*apply_slant + L1BS.wd_corr.*apply_wd;
%     [~,win_delay_ref_index_2] = min(L1BS.shift(:));
%     L1BS.shift(:) = L1BS.shift(:)- L1BS.shift(win_delay_ref_index_2);


L1BS.shift_coarse = round(L1BS.shift);
L1BS.shift_fine = L1BS.shift-L1BS.shift_coarse;


%% increase the window size 
if(avoid_wrapps_cnf)
    L1BS.N_windows  = ceil(max(max(abs(L1BS.shift))/N_samples)+1);
    truncated_zeros = zeros(1,N_samples*zp_fact_range_cnf*(L1BS.N_windows-1));
else
    L1BS.N_windows  = 1;
    truncated_zeros =[];
end

% L1BS.beam_geo_corr          = zeros(1,max(L1BS.N_beams_stack),N_samples*L1BS.N_windows);
% 
%  if(strcmp(mode,'SIN'))
%     L1BS.beam_geo_corr_2 	= zeros(1,max(L1BS.N_beams_stack),N_samples*L1BS.N_windows);
%  end



%% 5. APPLYING THE ALIGNMENTS IN TIME 
    
%     wfm_geo_corr_aux = zeros (N_max_beams_stack_chd,N_samples);
%     if(strcmp(mode,'SIN'))
%         wfm_geo_corr_aux_2 = zeros (N_max_beams_stack_chd,N_samples);
%     end
    
    
    
      
%         wfm_geo_corr_aux(i_beam,:) = squeeze(L1BS.beams_surf(i_beam,:)).' .* exp(2i*pi_cst/N_samples*L1BS.shift(i_beam).*i_samples);
%         beam_zp_fft                     = fftshift(fft(wfm_geo_corr_aux(i_beam,:).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples);
%         beam_surf_aux1_fft_increased    = (cat(2,beam_zp_fft,truncated_zeros));
%         beam_surf_aux1_fft_aligned      = circshift(beam_surf_aux1_fft_increased,[0,L1BS.shift_coarse(i_beam)*zp_fact_range_cnf]);
%         L1BS.beam_geo_corr(i_beam,:) = beam_surf_aux1_fft_aligned;
        L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:) = L1BS.beams_surf(1:L1BS.N_beams_stack,:) .* exp(2i.*pi_cst./N_samples.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples);        
        if(strcmp(mode,'SIN')) && strcmp(processing_mode_cnf,'SIN')
%          	wfm_geo_corr_aux_2(i_beam,:)    = squeeze(L1BS.beams_surf_2(i_beam,:)).' .* exp(2i*pi_cst/N_samples*L1BS.shift(i_beam).*i_samples);
%             beam_zp_fft                     = fftshift(fft(wfm_geo_corr_aux_2(i_beam,:).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples);
%             beam_surf_aux1_fft_increased    = (cat(2,beam_zp_fft,truncated_zeros));
%             beam_surf_aux1_fft_aligned      = circshift(beam_surf_aux1_fft_increased,[0,L1BS.shift_coarse(i_beam)*zp_fact_range_cnf]);
%             L1BS.beam_geo_corr_2(i_beam,:) = beam_surf_aux1_fft_aligned;
             L1BS.beam_geo_corr_2(1:L1BS.N_beams_stack,:)    = L1BS.beams_surf_2(1:L1BS.N_beams_stack,:) .* exp(2i.*pi_cst./N_samples.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples);
        end     
%        if abs(L1BS.shift_coarse(i_beam))* zp_fact_range_cnf >= N_samples* zp_fact_range_cnf
%          L1BS.good_samples(i_beam,:)=0;
%        else

%     for i_beam = 1:L1BS.N_beams_stack
%         if(L1BS.shift_coarse(i_beam)>0)
%             start_sample=ceil(L1BS.shift_coarse(i_beam))*zp_fact_range_cnf+1;
%             final_sample=N_samples*zp_fact_range_cnf;
%         else
%             start_sample=1;
%             final_sample=(N_samples*zp_fact_range_cnf+floor((L1BS.shift_coarse(i_beam))*zp_fact_range_cnf));
%         end
%         
%         L1BS.good_samples(i_beam,start_sample:final_sample)=1;
%     end

    for i_beam = 1:L1BS.N_beams_stack
        if(L1BS.shift_coarse(i_beam)>0 && L1BS.shift_coarse(i_beam)<N_samples)
            start_sample=ceil(L1BS.shift_coarse(i_beam))*zp_fact_range_cnf+1;
            final_sample=N_samples*zp_fact_range_cnf;
            L1BS.good_samples(i_beam,1:start_sample-1)=0;
            L1BS.good_samples(i_beam,start_sample:final_sample)=1;
        elseif(L1BS.shift_coarse(i_beam)<=0 && L1BS.shift_coarse(i_beam)>-1*N_samples)
            start_sample=1;
            final_sample=(N_samples*zp_fact_range_cnf+floor((L1BS.shift_coarse(i_beam))*zp_fact_range_cnf));
            L1BS.good_samples(i_beam,start_sample:final_sample)=1;
            L1BS.good_samples(i_beam,final_sample+1:N_samples*zp_fact_range_cnf)=0;
        else
            L1BS.good_samples(i_beam,1:N_samples*zp_fact_range_cnf)=0;
        end
    end
    

% L1BS = rmfield(L1BS,'beams_surf');
% if(strcmp(mode,'SIN'))
%     L1BS = rmfield(L1BS,'beams_surf_2');
% end
end

    


