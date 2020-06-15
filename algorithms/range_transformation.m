%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the RANGE COMPRESSION as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: The purpose of the range compression is to perform range 
% transformation of the input stacks (with an FFT) and then generate the 
% power waveforms.
%
% Revisions 
% 1.0 
% 2.0 Range transformation done per record. New architecture
% 2.1 Remove the fftshift on the S3 case
% 2.2 Stack Coherence & Phase diff flag 
% 2.3 Corrected 
% 2.4 added cnf flag processing_mode_cnf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L1BS] = range_transformation (L1BS)

global N_samples zp_fact_range_cnf coherence_threshold_cnf
global mode mission processing_mode_cnf
global window_rg_cnf window_range
global stack_interferometry_flag_cnf

% L1BS.beams_rng_cmpr = zeros (max(L1BS.N_beams_stack), N_samples * zp_fact_range_cnf*L1BS.N_windows);
% if(strcmp(mode,'SIN'))
%     L1BS.beams_rng_cmpr_2 = zeros (max(L1BS.N_beams_stack), N_samples * zp_fact_range_cnf*L1BS.N_windows);
%     L1BS.phase_diff     = zeros (max(L1BS.N_beams_stack), N_samples *zp_fact_range_cnf*L1BS.N_windows);
%     L1BS.coherence      = zeros (max(L1BS.N_beams_stack), N_samples * zp_fact_range_cnf*L1BS.N_windows);
%     L1BS.coherence_mask= zeros (max(L1BS.N_beams_stack), N_samples * zp_fact_range_cnf*L1BS.N_windows);
% end
% beams_rng_cmpr_I = zeros (L1BS.N_total_surf_loc, max(L1BS.N_beams_stack), N_samples * zp_fact_range_cnf);
% beams_rng_cmpr_Q = zeros (L1BS.N_total_surf_loc, max(L1BS.N_beams_stack), N_samples * zp_fact_range_cnf);
% beams_rng_cmpr_aux = zeros (L1BS.N_total_surf_loc, max(L1BS.N_beams_stack), N_samples * zp_fact_range_cnf);
   
    

    %% 1.FFT
	switch mission
		case 'CR2'
			if(window_rg_cnf)
				% Hamming        
				beams_zp_fft    = fftshift(fft((squeeze(L1BS.beam_geo_corr(:,:)).*(ones(length(L1BS.beam_geo_corr(:,1)),1)*window_range)).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples* zp_fact_range_cnf);
				if(strcmp(mode,'SIN')) && strcmp(processing_mode_cnf,'SIN') 
					beams_zp_fft_2  = fftshift(fft((squeeze(L1BS.beam_geo_corr_2(:,:)).*(ones(length(L1BS.beam_geo_corr(:,1)),1)*window_range)).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples* zp_fact_range_cnf);
				end
			else
				beams_zp_fft    = fftshift(fft(squeeze(L1BS.beam_geo_corr(:,:)).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples* zp_fact_range_cnf);
				if(strcmp(mode,'SIN')) && strcmp(processing_mode_cnf,'SIN')
					beams_zp_fft_2  = fftshift(fft(squeeze(L1BS.beam_geo_corr_2(:,:)).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples* zp_fact_range_cnf);
				end
			end
		case {'S3_','S3A','S3B'} % don't do the shift
			if(window_rg_cnf)
				% Hamming        
				beams_zp_fft    = fftshift(fft((squeeze(L1BS.beam_geo_corr(:,:)).*(ones(length(L1BS.beam_geo_corr(:,1)),1)*window_range)).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples* zp_fact_range_cnf);
			else
				beams_zp_fft    = fftshift(fft(squeeze(L1BS.beam_geo_corr(:,:)).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples* zp_fact_range_cnf);
			end
	end
     
     
       
     
    
    	%% 2. AMPLITUDE to POWER
        
        L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,:)   = abs(beams_zp_fft(1:L1BS.N_beams_stack,:)).^2;
        L1BS.beams_rng_cmprIQ(1:L1BS.N_beams_stack,:) = (beams_zp_fft(1:L1BS.N_beams_stack,:));
        
        if(strcmp(mode,'SIN')) && strcmp(processing_mode_cnf,'SIN')
            
            L1BS.beams_rng_cmpr_2(1:L1BS.N_beams_stack,:) = abs(beams_zp_fft_2(1:L1BS.N_beams_stack,:)).^2;
            L1BS.beams_rng_cmprIQ_2(1:L1BS.N_beams_stack,:)= (beams_zp_fft_2(1:L1BS.N_beams_stack,:));
            if(stack_interferometry_flag_cnf == 1)
                % substacking performed, otherwise coherence will be always 1
                for i_beam=1:L1BS.N_beams_stack
                    beam_margin_l=2;
                    beam_margin_r=2;
                    N_substack = beam_margin_l+ beam_margin_r+1;
                    if(i_beam-beam_margin_l<1)
                        beam_margin_l=i_beam-1;
                        beam_margin_r = N_substack-beam_margin_l-1;
                        continue; % do not compute anything in the edges
                        
                    elseif((i_beam+beam_margin_r)>L1BS.N_beams_stack)
                        beam_margin_r=L1BS.N_beams_stack-i_beam;
                        beam_margin_l=N_substack-beam_margin_r-1;
                        continue;
                    end
                    beams = (i_beam-beam_margin_l):(i_beam+beam_margin_r);
                    aux_1 = L1BS.beams_rng_cmprIQ(beams,:);
                    aux_2 = L1BS.beams_rng_cmprIQ_2(beams,:);
                    aux_1(find(aux_1==0)) = nan;
                    aux_2(find(aux_2==0)) = nan;
                    stack_power(i_beam,:) = nanmean(L1BS.beams_rng_cmprIQ(beams,:).*conj(L1BS.beams_rng_cmprIQ_2(beams,:)));
                    a1_power(i_beam,:) = (nanmean(L1BS.beams_rng_cmpr(beams,:)));
                    a2_power(i_beam,:) = (nanmean(L1BS.beams_rng_cmpr_2(beams,:)));
                    if(L1BS.Gap_flag(i_beam))
                        L1BS.coherence(i_beam,:) = (abs(stack_power(i_beam,:)))./sqrt(a1_power(i_beam,:).*a2_power(i_beam,:));
                        L1BS.coherence_mask(i_beam,find(L1BS.coherence(i_beam,:)>coherence_threshold_cnf)) = 1; % to be checked

                    end
%                     phase1 = atan(imag(beams_zp_fft(1:L1BS.N_beams_stack,:))./real(beams_zp_fft(1:L1BS.N_beams_stack,:)))*180/pi;
%                     phase2 = atan(imag(beams_zp_fft_2(1:L1BS.N_beams_stack,:))./real(beams_zp_fft_2(1:L1BS.N_beams_stack,:)))*180/pi;
%                     L1BS.phase_diff(1:L1BS.N_beams_stack,:) = phase1 - phase2;
%                     stack_power = beams_zp_fftIQ(1:L1BS.N_beams_stack,:).*conj(beams_zp_fftIQ_2(1:L1BS.N_beams_stack,:));
%                     a1_power = squeeze(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,:));
%                     a2_power = squeeze(L1BS.beams_rng_cmpr_2(1:L1BS.N_beams_stack,:));
%                     L1BS.coherence(1:L1BS.N_beams_stack,:) = sqrt((abs(stack_power).^2)./(a1_power.*a2_power));
%                     L1BS.coherence_mask(1:L1BS.N_beams_stack,find(L1BS.coherence(1:L1BS.N_beams_stack,:)>coherence_threshold_cnf)) = 1; % to be checked
                end
            %Based on the conventional CryoSAt-2 SARIn processing:
            %coherencey and phase info is after multilooking
            
%             
            end
                
%         %% 3. FORCING WRAPPED SAMPLES TO 0 (applying the mask)
%         if abs(L1BS.shift_coarse(i_beam)) < N_samples
%             if L1BS.shift_coarse(i_beam) > 0
%                 limit1 = 1;
%                 limit2 = zp_fact_range_cnf * L1BS.shift_coarse(i_beam);
%             else
%                 limit1 = (N_samples - abs(L1BS.shift_coarse(i_beam)))*zp_fact_range_cnf;
%                 limit2 = N_samples*zp_fact_range_cnf;
%             end
%         else
%             limit1 = 1;
%             limit2 = N_samples*zp_fact_range_cnf;
%         end
%         L1BS.beams_rng_cmpr(i_beam,limit1:limit2) = 0;
    
    

end

