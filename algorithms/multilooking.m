%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the MULTI-LOOKING as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: The purpose of the multi-looking is to performing an average 
% of all the beams and create the L1B waveform.
% 
%
% ----------------------------------------------------------
% Author:    Roger Escol�  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  M�nica Roca   / isardSAT
% Last rev.: M�nica Roca   / isardSAT (11/09/2013)
% 
% Version
% 1.0 2015/01/01 From S6
% 1.1 Added stack mask vector else condition 512
% 1.2 Added apply_stack_mask_cnf
% 2.0 Multilooking per record. New architecture. N_max_beams_stack added as
% global
% 2.1 For loop removed when averaging
% 2.2 Statistics moved after mask application
% 2.3 Changed N_samples_sar_chd by N_samples
% 2.4 L1BS.N_windows included in the initialisation of the L1B vectors
% 2.5 added cnf flag processing_mode_cnf
%
% The purpose of the multi-looking is to correct the beams for the window
% delay misalignments, creating the stack for every surface location and
% performing an average of all the waveforms in it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L1BS,L1B]   = multilooking(L1BS)

global N_samples pi_cst
global zp_fact_range_cnf N_beams_sub_stack_cnf use_zeros_cnf delete_ambiguities_cnf
global i_sample_start_chd compute_L1_stadistics_cnf coherency_mask_cnf apply_stack_mask_cnf
global mode N_max_beams_stack_chd processing_mode_cnf
global avoid_beams_mask_allzeros_cnf
global mask_look_angles_CR2_cnf mask_look_angles_CR2_min_chd mask_look_angles_CR2_max_chd
% global k_noise_top_cnf k_noise_floor_cnf
% global start_angle_spw weights_spw angles_spw range_index1_cnf range_index2_cnf
global avoid_noisy_beams_cnf method_noise_est_cnf method_noisy_thresholding_cnf noise_estimation_window param_std_noise_thres_cnf

N_sub_stacks = floor(L1BS.N_beams_stack / N_beams_sub_stack_cnf);
a_gauss                             = 0;
L1B.stack_look_ang_centre           = 0;
L1B.stack_std                       = 0;
L1B.stack_width                     = 0;
L1B.stack_pointing_ang_centre       = 0;
L1B.stack_skewness                  = 0;
L1B.stack_kurtosis                  = 0;
L1B.wfm_cor_i2q2_sar_ku             = zeros (1,N_samples*zp_fact_range_cnf*L1BS.N_windows);
if(strcmp(mode,'SIN')) && strcmp(processing_mode_cnf,'SIN')
    L1B.phase_diff                      = zeros (1,N_samples*zp_fact_range_cnf*L1BS.N_windows);
    L1B.coherence                       = zeros (1,N_samples*zp_fact_range_cnf*L1BS.N_windows);
end
wrap_margin         = 15; % range of samples to avoid considering part of the mask vector
rmc_margin          = 6;
rmc_mask            = ones  (N_max_beams_stack_chd,N_samples*zp_fact_range_cnf);
rmc_mask(:,(N_samples/2+i_sample_start_chd-1-rmc_margin)*zp_fact_range_cnf:N_samples*zp_fact_range_cnf)    = 0;
if i_sample_start_chd > 1
  rmc_mask(:,1:i_sample_start_chd-1)= 0;
end

L1B.stack_mask_vector           = ones (1,N_max_beams_stack_chd);%zeros (1,N_max_beams_stack_chd)+N_samples*zp_fact_range_cnf;
power_fitted                    = zeros (1,N_max_beams_stack_chd);
L1B.N_beams_used                = zeros (1,N_samples*zp_fact_range_cnf*L1BS.N_windows);
beams_off                       = 0;
doppler_central                 = 0;
stack_beam_centre               = 0;
L1B.N_beams_contributing        = 0;


    


%GAP mask
% for i_beam = 1:L1BS.N_beams_stack
%     L1BS.good_samples(i_beam,:)=L1BS.good_samples(i_beam,:) *L1BS.Gap_flag(i_beam);
% end
L1BS.good_samples=L1BS.good_samples.*(L1BS.Gap_flag.'*ones(1,N_samples*zp_fact_range_cnf));
%RMC mask
if L1BS.ProcessID_surf== 59;
  L1BS.good_samples(:,:) = squeeze(L1BS.good_samples(:,:)).*rmc_mask;
end

%CAL 4 Mask
cal_beams=find(L1BS.ProcessID_beam(:)== 57);
L1BS.good_samples(cal_beams,:)  = 0;

%% BUILD THE AMBIGUITY MASK


% diff_amb_range = diff(range_sat_surf_amb(100,:))* 2 / c_cst *zp_fact_range_cnf/ mean(T0_sar_surf(100,:));
% p=polyfit(1:length(diff_amb_range),diff_amb_range,1);
% 
% ambiguity_slope = abs(p(2)*1.35); % [samples/beam]   

% the 2 is because we have substracted the slant range unproperly for the ambiguous samples.
% before the geometry corrections the slope is just p(2) not 1.35*p(2). 

% The ambiguity slope should be aprox 3.66 (visualy computed)


if(delete_ambiguities_cnf==1)
    ambiguity_slope = 3.66;
    ambiguity_margin_back = 10*zp_fact_range_cnf; % Cut 8*zp samples before the ambiguity
    ambiguity_margin_forw = 30*zp_fact_range_cnf; % Cut 8*zp samples before the ambiguity
    stack_mispointing = 16; % Number of beams of mispointing (as the pointing beam index is the 33 at burst levels)

    surf_completed = find(abs(diff(diff(L1BS.N_beams_stack)))>2);
    % surf completed == 1 when the N_beams oscilates from +1 -1 of the N_max_beams_stack_chd
    % surf completed > 2 in the beginning and end of the product



    ambiguity_mask                  = ones  (N_max_beams_stack_chd,N_samples*zp_fact_range_cnf);

    % compute the leading edge position, cosidering only the 20 central beams
    [~,doppler_central]= min(abs(pi_cst/2-L1BS.beam_ang_surf(:)));
    if(doppler_central>L1BS.N_beams_stack-9)
        right_beams_cnf=L1BS.N_beams_stack-doppler_central;
    else    
        right_beams_cnf=9;
    end
    if(doppler_central<11)
        left_beams_cnf=doppler_central-1;
    else
        left_beams_cnf=10;
    end
    [~,max_pos_beams]=(max(squeeze(L1BS.beams_rng_cmpr(:,:)).'));
    max_pos = round(min(max_pos_beams(doppler_central-left_beams_cnf:doppler_central+right_beams_cnf)));


    % compute the offset to the ambiguity mask for the incomplete stacks at the end of the product
    beams_off = N_max_beams_stack_chd-L1BS.N_beams_stack;

    if(L1BS.surf_counter<surf_completed(1)) 
        beams_off_back=beams_off;
        beams_off_forw=0;
    elseif(i_surf>surf_completed(2))
        beams_off_forw=beams_off;
        beams_off_back=0;
    else
        beams_off_forw=0;
        beams_off_back=0;
    end

    %forward beams
    for i_beam = L1BS.N_beams_stack:-1:1 
        sample_amb_forw(i_beam) = max_pos - (-beams_off_forw-L1BS.N_beams_stack-stack_mispointing-1+i_beam)*ambiguity_slope -ambiguity_margin_forw;
        if(round(sample_amb_forw(i_beam))<N_samples*zp_fact_range_cnf)
            ambiguity_mask(i_beam,round(sample_amb_forw(i_beam)):end)=0;
        end
        if(round(sample_amb_forw(i_beam))<max_pos)
            ambiguity_mask(i_beam,:)=0;
        end
    end

    % not completed stacks
    for i_beam = L1BS.N_beams_stack+1:N_max_beams_stack_chd
        ambiguity_mask(i_beam,:) = 0;
    end

    %backward beams
    sample_amb_backw = zeros(1,N_max_beams_stack_chd);

    for i_beam = 1:L1BS.N_beams_stack
        sample_amb_backw(i_beam) = max_pos - (-beams_off_back+ stack_mispointing-i_beam) * ambiguity_slope-ambiguity_margin_back;
        if ((round(sample_amb_backw(i_beam))<N_samples*zp_fact_range_cnf) && round(sample_amb_backw(i_beam))>0)
            ambiguity_mask(i_beam,round(sample_amb_backw(i_beam)):end)=0;
        end
        if(round(sample_amb_backw(i_beam))<max_pos)
            ambiguity_mask(i_beam,:)=0;
        end
        L1BS.stack_mask(i_beam,:) = ambiguity_mask(i_beam,:).*squeeze(L1BS.good_samples(i_beam,:))';
    end

else

%     for i_beam = 1:L1BS.N_beams_stack
%         L1BS.stack_mask(i_beam,:) = L1BS.good_samples(i_beam,:);
%     end
      L1BS.stack_mask = L1BS.good_samples;
end

%include the additional mask depending on the look angles as in CR2
if(mask_look_angles_CR2_cnf)
   indexes_mask_look_angles=find(L1BS.look_ang_surf<=mask_look_angles_CR2_min_chd | L1BS.look_ang_surf>=mask_look_angles_CR2_max_chd);
   if ~isempty(indexes_mask_look_angles)
        L1BS.stack_mask(indexes_mask_look_angles,:)=0;
   end
end

% ----------------- mask out those outer noisy beams ----------------------
if avoid_noisy_beams_cnf
    %In a first approach not filtering any beam depending on the mask
    valid_noise_samples=noise_estimation_window(1)*zp_fact_range_cnf:noise_estimation_window(2)*zp_fact_range_cnf;
    switch method_noise_est_cnf
        case 'before_geom_corr'
            beams_rng_cmpr_before_geom_corr=abs(fftshift(fft(squeeze(L1BS.beams_surf(:,:)).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples* zp_fact_range_cnf)).^2;
            if(use_zeros_cnf ==0)
                valid_samples=reshape(beams_rng_cmpr_before_geom_corr(1:L1BS.N_beams_stack,valid_noise_samples),[1,L1BS.N_beams_stack*length(valid_noise_samples)]);
            elseif(use_zeros_cnf==1)
                valid_samples=reshape(beams_rng_cmpr_before_geom_corr(1:L1BS.N_beams_stack,valid_noise_samples).*L1BS.stack_mask(1:L1BS.N_beams_stack,valid_noise_samples),...
                    [1,L1BS.N_beams_stack*length(valid_noise_samples)]);
            end
        case 'after_geom_corr'
            if(use_zeros_cnf ==0)
                valid_samples=reshape(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,valid_noise_samples),[1,L1BS.N_beams_stack*length(valid_noise_samples)]);
            elseif(use_zeros_cnf==1)
                valid_samples=reshape(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,valid_noise_samples).*L1BS.stack_mask(1:L1BS.N_beams_stack,valid_noise_samples),...
                    [1,L1BS.N_beams_stack*length(valid_noise_samples)]);
            end
    end
    
    mean_noise=mean(valid_samples); 
    std_noise=std(valid_samples,1); %std as = 1/N*(sum_i=0^N-1(x_i-mean_x)^2)
    noisy_wvfms=mean(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,valid_noise_samples),2)>(mean_noise+param_std_noise_thres_cnf*std_noise);
    if any(noisy_wvfms)
        L1BS.stack_mask(noisy_wvfms,:)=0;
    end
    
    %should try to use the same set of samples for all beams
%         valid_noise_samples=zeros(L1BS.N_beams_stack,N_samples*zp_fact_range_cnf);
%     valid_noise_samples(:,noise_estimation_window(1)*zp_fact_range_cnf:noise_estimation_window(2)*zp_fact_range_cnf)=1;
%     valid_noise_samples=find(prod(valid_noise_samples,1));
%    valid_noise_samples=find(prod(valid_noise_samples.*L1BS.stack_mask(1:L1BS.N_beams_stack,:),1));        
    % noise statistics
%     switch method_noise_est_cnf
%         case 'Beam_based'
%             %Noise estimated per beam or look
%             mean_noise=mean(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,valid_noise_samples),2);
%             std_noise=std(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,valid_noise_samples),0,2);
%         case 'Stack_based'
%             %Noise estimated using the whole stack information
%             valid_samples=reshape(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,valid_noise_samples),[1,L1BS.N_beams_stack*length(valid_noise_samples)]);
%             mean_noise=mean(valid_samples,1)*ones(L1BS.N_beams_stack,1);
%             std_noise=std(valid_samples,0)*ones(L1BS.N_beams_stack,1);
%     end
%     %avoid those samples set to 0 --> NaN
%     beams_rng_cmpr=L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,:);
%     beams_rng_cmpr(~L1BS.stack_mask(1:L1BS.N_beams_stack,:))=nan;
%     switch method_noisy_thresholding_cnf
%         case 'Average_based'
%             %Average value within waveform for comparison against
%             %threshold
%             wvfm_value_thres=nanmean(beams_rng_cmpr,2);            
%         case 'Peak_based'
%             %Peak value for comparison against threshold
%             wvfm_value_thres=nanmax(beams_rng_cmpr,[],2);
%     end
%     noisy_wvfms=wvfm_value_thres<=(mean_noise+param_std_noise_thres_cnf*std_noise);
%     %remove/mask those noisy beams
%     if any(noisy_wvfms)
%         L1BS.stack_mask(noisy_wvfms,:)=0;
%     end
    
end

if(apply_stack_mask_cnf)

%     for i_beam = 1:L1BS.N_beams_stack
%         find_zeros = find((L1BS.stack_mask(i_beam,:))==0);
%         if(~isempty(find_zeros))
%             find_jumps = find(diff(find_zeros)>1); 
%             % findjumps should be only one value unless strage maskings are applied in the future
%             if(isempty(find_jumps))
%                 find_jumps=0;
%             end
%             aux=min(find_zeros(find(find_zeros==find_zeros(find_jumps+1)):end));
%             if length(find_zeros)<wrap_margin && aux<wrap_margin
%                 L1B.stack_mask_vector (i_beam) = N_samples*zp_fact_range_cnf*L1BS.N_windows;
%             else
%                 L1B.stack_mask_vector (i_beam) = aux;
%             end
%         else
%             L1B.stack_mask_vector (i_beam) = N_samples*zp_fact_range_cnf*L1BS.N_windows; 
%         end
%     end
    for i_beam = 1:L1BS.N_beams_stack
        if L1BS.stack_mask(i_beam,end)~=1
            a=find(L1BS.stack_mask(i_beam,:),1,'last');
            if ~isempty(a)
                L1B.stack_mask_vector (i_beam)=a+1;
            else
                L1B.stack_mask_vector (i_beam)=1;
            end
        else
            L1B.stack_mask_vector (i_beam)=N_samples.*zp_fact_range_cnf;
        end
    end
else
    L1B.stack_mask_vector (1:L1BS.N_beams_stack) = N_samples*zp_fact_range_cnf;
end

% a bit tricky. Divide or Multiply by the mask in order to force the
% aliased samples to infinity or to 0 so then are considered or not in
% the mean operation. Ask Albert for more explanation. 
if avoid_beams_mask_allzeros_cnf
    %mask: beams where all the samples are set to 0 are non-contributing to the
    %ML
    idx_rows_int=any(L1BS.stack_mask,2);
    if any(idx_rows_int)
        L1B.start_beam = find(idx_rows_int,1,'first'); %first beam where not all samples are set to zero
        L1B.stop_beam = find(idx_rows_int,1,'last'); %last beam where not all samples are set to zero        
    else
        %all beams masked out no contribution is expected
        L1B.start_beam = 0; %first beam where not all samples are set to zero
        L1B.stop_beam = 1; %last beam where not all samples are set to zero       
    end
    
    L1B.N_beams_contributing = L1BS.N_beams_stack-length(find(L1B.stack_mask_vector(1:L1BS.N_beams_stack)==1));
else
    idx_rows_int=1:1:L1BS.N_beams_stack;
    L1B.start_beam=1;
    L1B.stop_beam=L1BS.N_beams_stack;
    L1B.N_beams_contributing = L1BS.N_beams_stack;
end
L1B.N_beams_start_stop = L1B.stop_beam-L1B.start_beam+1; % number of beams or looks from start to stop
%move all the contributing beams to the very beginning of stack mask vector
% using a similar approach as proposed for Sentinel-6
if L1B.start_beam~=0
    aux=L1B.stack_mask_vector(L1B.start_beam:L1B.stop_beam);
    L1B.stack_mask_vector(1:1:L1BS.N_beams_stack)=1;
    L1B.stack_mask_vector(1:L1B.N_beams_start_stop)=aux;
else
    L1B.stack_mask_vector(1:1:L1BS.N_beams_stack)=1;
end

if(apply_stack_mask_cnf)
    if(use_zeros_cnf ==0)
        L1BS.beams_rng_cmpr (:,:) = L1BS.beams_rng_cmpr(:,:)./L1BS.stack_mask(:,:);
        L1BS.beams_rng_cmprIQ (:,:) = L1BS.beams_rng_cmprIQ(:,:)./L1BS.stack_mask(:,:);
    elseif(use_zeros_cnf==1)
        L1BS.beams_rng_cmpr (:,:) = L1BS.beams_rng_cmpr(:,:).*L1BS.stack_mask(:,:);
        L1BS.beams_rng_cmprIQ (:,:) = L1BS.beams_rng_cmprIQ(:,:).*L1BS.stack_mask(:,:);
    end      
end

%% 3. STACK PARAMETERS (need to be updated considering only rows with not all elements to zero)
%Power of the beams
if(compute_L1_stadistics_cnf)
    central_beams = 0;
    beam_power = zeros (1,L1BS.N_beams_stack);
    for i_beam = 1:L1BS.N_beams_stack
        beam_power(i_beam) = sum (L1BS.beams_rng_cmpr (i_beam,:)); 
    end
    beam_power=beam_power/max(beam_power);

    %Sub-Stacking
    power_sub_stack = zeros(N_sub_stacks);
    for i_beam = 1:N_sub_stacks
        power_sub_stack(i_beam) = sum (beam_power(((i_beam-1)*N_beams_sub_stack_cnf)+1:((i_beam-1)*N_beams_sub_stack_cnf)+N_beams_sub_stack_cnf)) / N_beams_sub_stack_cnf;
    end
    [~,doppler_central]= min(abs(pi_cst/2-L1BS.beam_ang_surf(:)));
    right_beams_cnf=124;
    if(doppler_central<125)
        left_beams_cnf=doppler_central-1;
    else
        left_beams_cnf=124;
    end

    if((doppler_central+right_beams_cnf)>L1BS.N_beams_stack)
        right_beams_cnf=L1BS.N_beams_stack-doppler_central;
    end
    central_beams_index = (doppler_central-left_beams_cnf:doppler_central+right_beams_cnf);
    valid_beams     = find(L1BS.Gap_flag(doppler_central-left_beams_cnf:doppler_central+right_beams_cnf));
    central_beams = central_beams_index(valid_beams);
    try
    % Gaussian fittting
        cfun_look                       = fit (L1BS.look_ang_surf(central_beams)', beam_power(central_beams)', 'gauss1'); %fitting vs look_angle
        cfun_pointing                   = fit (L1BS.pointing_ang_surf(central_beams)', beam_power(central_beams)', 'gauss1'); %fitting vs pointing_angle
        a_gauss                 = cfun_look.a1;
        L1B.stack_look_ang_centre       = cfun_look.b1;
        L1B.stack_pointing_ang_centre   = cfun_pointing.b1;
        L1B.stack_std               = cfun_look.c1/2;
        L1B.stack_width             = 2*sqrt(2*log(2))*cfun_look.c1;
        power_fitted(central_beams) = a_gauss * exp (-(L1BS.look_ang_surf(central_beams) - L1B.stack_look_ang_centre).^2 /(2*L1B.stack_std.^2));
        % Compute the characterization parameters:
        L1B.stack_skewness= skewness(power_fitted(central_beams));
        L1B.stack_kurtosis= kurtosis(power_fitted(central_beams))-3;

        cfun_beams = fit ((central_beams)', beam_power(central_beams)', 'gauss1'); %fitting vs look_angle
        L1B.stack_beam_centre       = cfun_beams.b1;
    catch
        L1B.stack_look_ang_centre       = 0;
        L1B.stack_pointing_ang_centre   = 0;
        L1B.stack_std               = 0;
        L1B.stack_width             = 0;
        L1B.stack_skewness= 0;
        L1B.stack_kurtosis= 0;
        L1B.stack_beam_centre       = 0;

    end
end


if(strcmp(mode,'SIN'))&& strcmp(processing_mode_cnf,'SIN')
    if(apply_stack_mask_cnf)
        if(use_zeros_cnf ==0)

            L1BS.beams_rng_cmpr_2 (:,:) = L1BS.beams_rng_cmpr_2(:,:)./L1BS.stack_mask(:,:);
            L1BS.beams_rng_cmprIQ_2 (:,:) = L1BS.beams_rng_cmprIQ_2(:,:)./L1BS.stack_mask(:,:);
% !!!! TO BE REVIEWED
%             L1BS.phase_diff (:,:) = L1BS.phase_diff(:,:)./L1BS.stack_mask(:,:);
%             L1BS.coherence (:,:) = L1BS.coherence(:,:)./L1BS.stack_mask(:,:);

%             if (coherency_mask_cnf)
%                 L1BS.beams_rng_cmpr_Coh (:,:) = L1BS.beams_rng_cmpr(:,:)./L1BS.coherence_mask(:,:);
%                 L1BS.beams_rng_cmpr_2_Coh (:,:) = L1BS.beams_rng_cmpr_2(:,:)./L1BS.coherence_mask(:,:);
%                 L1BS.beams_rng_cmprIQ_Coh  (:,:) = L1BS.beams_rng_cmprIQ  (:,:)./L1BS.coherence_mask(:,:);
%                 L1BS.beams_rng_cmprIQ_2_Coh(:,:) = L1BS.beams_rng_cmprIQ_2(:,:)./L1BS.coherence_mask(:,:);
%                 L1BS.phase_diff_Coh (:,:)        = L1BS.phase_diff(:,:)./L1BS.coherence_mask(:,:);
%             end
        elseif(use_zeros_cnf==1)
            L1BS.beams_rng_cmpr_2 (:,:) = L1BS.beams_rng_cmpr_2(:,:).*L1BS.stack_mask(:,:);
%             L1BS.phase_diff (:,:) = L1BS.phase_diff(:,:).*L1BS.stack_mask(:,:);
%             L1BS.coherence (:,:) = L1BS.coherence(:,:).*L1BS.stack_mask(:,:);
        end
    end

end

%multilooking /averaging 
L1BS.beams_rng_cmpr         (~isfinite(L1BS.beams_rng_cmpr))        = NaN;
L1B.N_beams_used            (1:N_samples*zp_fact_range_cnf) = (sum(isfinite(L1BS.beams_rng_cmpr(idx_rows_int,:))));
L1B.wfm_cor_i2q2_sar_ku     (1:N_samples*zp_fact_range_cnf) = nanmean(L1BS.beams_rng_cmpr(idx_rows_int,:),1);
L1BS.beams_rng_cmprIQ (~isfinite(L1BS.beams_rng_cmprIQ))      = NaN;


if(strcmp(mode,'SIN')) && strcmp(processing_mode_cnf,'SIN')
    L1BS.beams_rng_cmpr_2       (~isfinite(L1BS.beams_rng_cmpr_2))      = NaN;
    L1B.wfm_cor_i2q2_sar_ku_2   (1:N_samples*zp_fact_range_cnf*L1BS.N_windows) = nanmean(L1BS.beams_rng_cmpr_2 (idx_rows_int,:),1);
    L1BS.beams_rng_cmprIQ_2 (~isfinite(L1BS.beams_rng_cmprIQ))      = NaN;
    %compute the L1B phase and coherence as conventional way using
    %multilook over the crossproduct and as noted by Wingham in 2006 paper
    cross_product=L1BS.beams_rng_cmprIQ.*conj(L1BS.beams_rng_cmprIQ_2);
    cross_product=nanmean(cross_product(idx_rows_int,:),1);
    L1B.phase_difference =angle(cross_product); %recoverying phase between [-pi,pi]
    L1B.coherence=abs(cross_product)./sqrt(L1B.wfm_cor_i2q2_sar_ku.*L1B.wfm_cor_i2q2_sar_ku_2);
    %due to possible masking some samples may go to zero after multilooking
    %this will lead to infinity values --> force them to zero if NaN keep
    %them
    L1B.coherence(isinf(L1B.coherence))=0;
    
    
    
% !!!!!! TO BE REVIEWED     
%     L1B.phase_diff_2          	(1:N_samples*zp_fact_range_cnf*L1BS.N_windows) = nanmean(L1BS.phase_diff       (idx_rows_int,:),1);
%     L1B.coherence_2             (1:N_samples*zp_fact_range_cnf*L1BS.N_windows) = nanmean(L1BS.coherence        (idx_rows_int,:),1);
% 
%     if (coherency_mask_cnf)
%         L1B.wfm_cor_i2q2_sar_kuIQ_Coh  (i_sample)    = mean(L1BS.beams_rng_cmprIQ_Coh(isfinite(L1BS.beams_rng_cmpr_Coh(idx_rows_int,i_sample)),i_sample));
%         L1B.wfm_cor_i2q2_sar_kuIQ_2_Coh  (i_sample)  = mean(L1BS.beams_rng_cmprIQ_2_Coh(isfinite(L1BS.beams_rng_cmpr_Coh(idx_rows_int,i_sample)),i_sample));
%         L1B.wfm_cor_i2q2_sar_ku_Coh       (i_sample) = mean(L1BS.beams_rng_cmpr_Coh(isfinite(L1BS.beams_rng_cmpr_Coh(idx_rows_int,i_sample)),i_sample));
%         L1B.wfm_cor_i2q2_sar_ku_2_Coh(i_sample)      = mean(L1BS.beams_rng_cmpr_2_Coh(isfinite(L1BS.beams_rng_cmpr_2_Coh(idx_rows_int,i_sample)),i_sample));
%     end
end

% !!!!!! TO BE REVIEWED
% if(strcmp(mode,'SIN'))
%     L1B.wfm_cor_i2q2_sar_kuIQ   (1:N_samples*zp_fact_range_cnf*L1BS.N_windows) = nanmean(L1BS.beams_rng_cmprIQ(idx_rows_int,:),1);
%     L1B.wfm_cor_i2q2_sar_kuIQ_2 (1:N_samples*zp_fact_range_cnf*L1BS.N_windows) = nanmean(L1BS.beams_rng_cmprIQ_2(idx_rows_int,:),1);
% 
% 
%     phase1 = atan(imag(L1B.wfm_cor_i2q2_sar_kuIQ(:))./real(L1B.wfm_cor_i2q2_sar_kuIQ(:)))*180/pi;
%     phase2 = atan(imag(L1B.wfm_cor_i2q2_sar_kuIQ_2(:))./real(L1B.wfm_cor_i2q2_sar_kuIQ_2(:)))*180/pi;
%     L1B.phase_diff(:) = (phase1 - phase2)';
% 
%     stack_power = L1B.wfm_cor_i2q2_sar_kuIQ (:).*conj(L1B.wfm_cor_i2q2_sar_kuIQ_2 (:)).^2;
%     a1_power    = L1B.wfm_cor_i2q2_sar_ku(:);
%     a2_power    = L1B.wfm_cor_i2q2_sar_ku_2(:);
%     L1B.coherence (:) = sqrt((abs(stack_power).^2)./(a1_power.*a2_power));
%     if(coherency_mask_cnf)
%         stack_power = L1B.wfm_cor_i2q2_sar_kuIQ_Coh (:).*conj(L1B.wfm_cor_i2q2_sar_kuIQ_2_Coh (:)).^2;
%         a1_power    = L1B.wfm_cor_i2q2_sar_ku_Coh(:);
%         a2_power    = L1B.wfm_cor_i2q2_sar_ku_2_Coh(:);
%         L1B.coherence_Coh (:) = sqrt((abs(stack_power).^2)./(a1_power.*a2_power));
%     end
% end
end