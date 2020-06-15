%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop
% This code implements the ACDC at stack level based on the GRSL "Amplitude
% and Dilation Compensation of the SAR Altimeter Backscattered Power"
%
% ---------------------------------------------------------
% Objective: Provide a flattened version of the stack a la ACDC and the
% corresponding multilooked ACDC waveform
% 
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Roger Escolà / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (02/06/2015)
% 
% 
% Last revision:    Eduard Makhoul / isardSAT V1 17/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -L1BS    =  structure containing the current stack information
%       -L1B     = 	structure containing the current L1B waveform
%       information
%       -
%      OPTIONAL
%       -LUT_f0_file =   file containing the look up table for the f0
%                       function parametrization (full path)
%       -LUT_f1_file =   file containing the look up table for the f1
%                       function parametrization (full path)
%       -path_Results=  full path to the folder results
%       
% OUTPUT:
%       -L1BS        =   structure of data with information at
%       stack level (include the ACDC stack)
%       -L1B         =   structure of data with information at level 1B (include the ACDC waveform)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - 
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0:
% V1.1: 
function [L1B,geo_param_buffer,fit_params_ini,accumulated,accumulated_error_epoch,accumulated_SSH,i_surf_fitted,flag_exit]  = ACDC_stack (L1A_buffer,L1BS,L1B,...
                                                                                                   range_index,delta_range,i_surf_stacked,...
                                                                                                   geo_param_buffer,fit_params_ini,...
                                                                                                   accumulated,accumulated_error_epoch,accumulated_SSH,i_surf_fitted,...
                                                                                                   func_f0,func_f1,files)
                                                                                               
disp(strcat('ACDC on surf #',{' '},num2str(i_surf_stacked)));                                                                                                                                                                                       
global cnf_p_ACDC
global apply_stack_mask_cnf
global wv_length_ku use_zeros_cnf
global N_ku_pulses_burst_chd N_samples zp_fact_range_cnf c_cst
global mission file_ext_string optional_ext_file_flag

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if nargin < 12
    func_f1=0; % initialize altitude to 0;
end
if nargin < 11
    func_f0=0; % initialize altitude to 0;
end

%% ----------------- CHECKING L1B WAVEFORM --------------------------------
if ~any(~isnan(L1B.wfm_cor_i2q2_sar_ku))
    flag_exit=-1;
    L1B.ACDC.waveform=NaN(1,N_samples*zp_fact_range_cnf);
    L1B.ACDC.range_index=NaN(1,N_samples*zp_fact_range_cnf);
    L1B.ACDC.Hs=NaN;
    L1B.ACDC.epoch=NaN;
    L1B.ACDC.Pu=NaN;
    return;
else
    flag_exit=0;
end

%% ----------------- INITIZALIZATION OF VARIABLES -------------------------
%--------------------------------------------------------------------------
%load in a structure of non fitting parameters
nf_p=gen_nonfit_params_ACDC(L1BS,L1B);

% ---------------- Stack + mask -------------------------------------------
%take only those ones within the start and stop
if(apply_stack_mask_cnf)
    stack_mask=L1BS.stack_mask(L1B.start_beam:L1B.stop_beam,:);
    if(use_zeros_cnf ==0)
        stack_mask(stack_mask==0)=NaN;        
    end
else
    stack_mask(:,:)=ones(nf_p.Neff,N_samples*zp_fact_range_cnf);
end
stack=L1BS.beams_rng_cmpr(L1B.start_beam:L1B.stop_beam,:).*stack_mask;
%for verbose purposes compute the stack before SRC
stack_before_SRC  = (abs(fftshift(fft(squeeze(L1BS.beams_surf(L1B.start_beam:L1B.stop_beam,:)).', N_samples * zp_fact_range_cnf),1).'/sqrt(N_samples* zp_fact_range_cnf))).^2.*stack_mask;

%------------------------- Look indexation --------------------------------
delta_look_angle=asin(wv_length_ku./(nf_p.pri_surf.*2.0*N_ku_pulses_burst_chd*nf_p.v_sat)); 
looks=(L1BS.look_ang_surf(L1B.start_beam:L1B.stop_beam)./delta_look_angle).'; %looks indexes   

%------------------ Method-dependent Initialization ---------------------
%--------------------------------------------------------------------------
switch cnf_p_ACDC.method_pre_fitting        
    case 0
        start_beam=L1B.start_beam;
        stop_beam=L1B.stop_beam;
    case 1
        % Specular beam (look for the beam closest to zero-Doppler)
        [~,beam_specular_point]=min(L1BS.doppler_ang_surf(L1B.start_beam:L1B.stop_beam));
        start_beam=max(beam_specular_point(1)-cnf_p_ACDC.num_left_beams_est_cnf+1,1);
        stop_beam=min(beam_specular_point(1)-cnf_p_ACDC.num_right_beams_est_cnf,L1B.N_beams_start_stop);
        nf_p.Neff=stop_beam-start_beam+1;
        stack=stack(start_beam:stop_beam,:);
        stack_mask=stack_mask(start_beam:stop_beam,:);
        stack_before_SRC=stack_before_SRC(start_beam:stop_beam,:);
        looks=looks(start_beam:stop_beam);
    case 2
        averaged_pow_act=nanmean(stack,2);
        [~,beam_max_power]=nanmax(averaged_pow_act);        
        start_beam=max(beam_max_power(1)-cnf_p_ACDC.num_left_beams_est_cnf+1,1);
        stop_beam=min(beam_max_power(1)-cnf_p_ACDC.num_right_beams_est_cnf,L1B.N_beams_start_stop);
        nf_p.Neff=stop_beam-start_beam+1;
        stack=stack(start_beam:stop_beam,:);
        stack_mask=stack_mask(start_beam:stop_beam,:);
        stack_before_SRC=stack_before_SRC(start_beam:stop_beam,:);
        looks=looks(start_beam:stop_beam);
end

%% ---------------- First estimation of Hs & Koff(epoch) ------------------
%--------------------------------------------------------------------------
%Preliminary estimation 
%if i_surf_stacked==1 || (mod(i_surf_stacked,cnf_p_ACDC.preliminary_est_surf_downsampling)==0)
    [Hs, epoch, a,corr_coef, ~, nf_p, flag,ml_wav_preliminary,ml_wav_fitted_preliminary]  = stack_fitting(stack,stack_mask,looks,range_index,...
                                                                              nf_p,cnf_p_ACDC,fit_params_ini(1:3),i_surf_stacked,...
                                                                              func_f0,func_f1);
    if flag==-1
        flag_exit=flag;
        L1B.ACDC.waveform=NaN(1,N_samples*zp_fact_range_cnf);
        L1B.ACDC.ml_wav_fitted_ACDC=NaN(1,N_samples*zp_fact_range_cnf);
        L1B.ACDC.ml_wav_preliminary=NaN(1,N_samples*zp_fact_range_cnf);
        L1B.ACDC.ml_wav_fitted_preliminary=NaN(1,N_samples*zp_fact_range_cnf);
        L1B.ACDC.range_index=NaN(1,N_samples*zp_fact_range_cnf);
        L1B.ACDC.Hs=NaN;
        L1B.ACDC.Hs_preliminary=NaN;
        L1B.ACDC.epoch=NaN;
        L1B.ACDC.epoch_preliminary=NaN;
        L1B.ACDC.SSH=NaN;
        L1B.ACDC.SSH_preliminary=NaN;
        L1B.ACDC.retracking_cor=NaN;
        L1B.ACDC.retracking_cor_preliminary=NaN;
        L1B.ACDC.tracker_range=NaN;
        L1B.ACDC.tracker_range_preliminary=NaN;
        L1B.ACDC.corr_coeff=NaN;
        L1B.ACDC.Pu=NaN;
        L1B.ACDC.sigma0=NaN;
        return
    else
        L1B.ACDC.ml_wav_preliminary=ml_wav_preliminary;
        L1B.ACDC.ml_wav_fitted_preliminary=ml_wav_fitted_preliminary;
        L1B.ACDC.Hs_preliminary=Hs;
        L1B.ACDC.epoch_preliminary=epoch;
        switch cnf_p_ACDC.range_index_method
            case 'conventional'
                L1B.ACDC.retracking_cor_preliminary  =   -(N_samples*zp_fact_range_cnf/2 - L1B.ACDC.epoch_preliminary)*delta_range;
            case 'resampling'
                L1B.ACDC.retracking_cor_preliminary  =   -(N_samples/2 - L1B.ACDC.epoch_preliminary)*delta_range;
        end
        L1B.ACDC.tracker_range_preliminary=L1BS.win_delay_surf*c_cst/2+L1B.ACDC.retracking_cor_preliminary;    
        %--------- Apply the geophysical corrections --------------------------
        if cnf_p_ACDC.geo_corr_application_flag
            [res] = geo_corrections_computation_per_surface(L1A_buffer,L1BS,cnf_p_ACDC);
            L1B.ACDC.SSH_preliminary=L1BS.alt_sat-(L1B.ACDC.tracker_range_preliminary+res.total_geo_corr);
        else
            L1B.ACDC.SSH_preliminary=L1BS.alt_sat-L1B.ACDC.tracker_range_preliminary;
        end
    end
%using the first estimate

if i_surf_stacked~=1
    %use the previous estimates to feed the ACDC
    %epoch=epoch+fit_params_ini(4);    
    Hs=fit_params_ini(2)*4;   
    
    %initialize the epoch using SSH previous estimate (within an sliding windowing)
    [res] = geo_corrections_computation_per_surface(L1A_buffer,L1BS,cnf_p_ACDC);
    retracker_corr=(L1BS.alt_sat-res.total_geo_corr-fit_params_ini(5)-L1BS.win_delay_surf*c_cst/2);
    switch cnf_p_ACDC.range_index_method
        case 'conventional'
            epoch  =1.0*retracker_corr/delta_range+N_samples*zp_fact_range_cnf/2;
        case 'resampling'
            epoch  =1.0*retracker_corr/delta_range+N_samples/2;
    end
    
end

for i_iter=1:cnf_p_ACDC.num_iterations
    %% --------------- Amplitude Compensation ---------------------------------
    %--------------------------------------------------------------------------
    [stack_ac,k,l,g]=amplitude_compensation(stack,looks,range_index,nf_p,cnf_p_ACDC,Hs,epoch);
    %k: includes already the epoch: range_index-epoch;
    
    stack_ac_before_smooth=stack_ac;
    
    %% -------------- SMOOTHING OF THE STACK ----------------------------------
    stack_ac_smooth=stack_ac;
    if cnf_p_ACDC.stack_smooth_active
        idx_rows_int=(any(stack_mask,2)); %avoid all-zeros or NaNs
        for i_beam=1:nf_p.Neff
            if idx_rows_int(i_beam)
                start=max(i_beam-cnf_p_ACDC.stack_smooth_win,1);
                final=min(i_beam+cnf_p_ACDC.stack_smooth_win,nf_p.Neff);
                stack_ac_smooth(i_beam,:)=nanmean(stack_ac(start:final,:),1);
            end
            
        end
        stack_ac_smooth=stack_ac_smooth.*stack_mask;
        stack_ac=stack_ac_smooth;
        %stack_ac=stack_ac_smooth;
    end
    
    %% ---------------- Dilation Compensation ---------------------------------
    %--------------------------------------------------------------------------
    %selection of the Doppler beam zero
    % norm_vel_sat=sqrt(L1BS.x_vel_sat_beam(start_beam:stop_beam).^2+ L1BS.y_vel_sat_beam(start_beam:stop_beam).^2+ L1BS.z_vel_sat_beam(start_beam:stop_beam).^2);
    % beam_angles=L1BS.beam_ang_surf(start_beam:stop_beam);
    % fd=-2/wv_length_ku.*norm_vel_sat.*cos(beam_angles);
    % [~,look_indexation_ref]=min(abs(fd)); %forcing it to 0.0 zero Doppler waveform
    %                          %This could be reviewed
    % look_indexation_ref=looks(look_indexation_ref(1));
    look_indexation_ref=0;
    
    [k_scaled]=dilation_compensation(k,g,nf_p,Hs,look_indexation_ref);
    %dilation_compensation(k,g,nf_p,SWH)
    %clear range_index;
    
    
    %% --------------- Retracking ACDC waveform -------------------------------
    %--------------------------------------------------------------------------
    %fit_params_ini(4)
    [waveform,k_ml_ACDC,estimates,flag,ml_wav_fitted_ACDC]=ACDC_retracking(stack_ac,k_scaled,nf_p,cnf_p_ACDC,[epoch,Hs,1.0,0.0],look_indexation_ref,i_surf_stacked,func_f0);
    if flag==-1
        L1B.ACDC.waveform=NaN(1,N_samples*zp_fact_range_cnf);
        L1B.ACDC.ml_wav_fitted_ACDC=NaN(1,N_samples*zp_fact_range_cnf);
        L1B.ACDC.range_index=NaN(1,N_samples*zp_fact_range_cnf);
        L1B.ACDC.Hs=NaN;
        L1B.ACDC.epoch=NaN;
        L1B.ACDC.SSH=NaN;
        L1B.ACDC.retracking_cor=NaN;
        L1B.ACDC.tracker_range=NaN;
        L1B.ACDC.corr_coeff=NaN;
        L1B.ACDC.Pu=NaN;
        L1B.ACDC.sigma0=NaN;
        
        %update the buffer: might be used for computation of initial params in a
        %sliding window approach
        geo_param_buffer(i_surf_stacked).epoch=L1B.ACDC.epoch;
        geo_param_buffer(i_surf_stacked).SSH=L1B.ACDC.SSH;
        if cnf_p_ACDC.rou_flag
            geo_param_buffer(i_surf_stacked).rou=estimates(2);
        else
            geo_param_buffer(i_surf_stacked).Hs=L1B.ACDC.Hs;
        end
        geo_param_buffer(i_surf_stacked).Pu=L1B.ACDC.Pu;
        geo_param_buffer(i_surf_stacked).error_epoch_ACDC=NaN;
        
        return
    else
        i_surf_fitted=i_surf_fitted+1;
        L1B.ACDC.waveform=waveform;
        L1B.ACDC.ml_wav_fitted_ACDC=ml_wav_fitted_ACDC;
        L1B.ACDC.range_index=k_ml_ACDC;
        L1B.ACDC.epoch=estimates(1);
        L1B.ACDC.Hs=estimates(2);
        L1B.ACDC.Pu=estimates(3);
        L1B.ACDC.sigma0=L1B.ACDC.Pu+L1B.wfm_scaling_factor_sar_ku;
        L1B.ACDC.corr_coeff=estimates(4);
        L1B.ACDC.flag_fitting=estimates(5);
        switch cnf_p_ACDC.range_index_method
            case 'conventional'
                L1B.ACDC.retracking_cor  =   -(N_samples*zp_fact_range_cnf/2 - L1B.ACDC.epoch)*delta_range;
            case 'resampling'
                L1B.ACDC.retracking_cor  =   -(N_samples/2 - L1B.ACDC.epoch)*delta_range;
        end
        L1B.ACDC.tracker_range=L1BS.win_delay_surf*c_cst/2+L1B.ACDC.retracking_cor;
        %--------- Apply the geophysical corrections --------------------------
        if cnf_p_ACDC.geo_corr_application_flag
            [res] = geo_corrections_computation_per_surface(L1A_buffer,L1BS,cnf_p_ACDC);
            L1B.ACDC.SSH=L1BS.alt_sat-(L1B.ACDC.tracker_range+res.total_geo_corr);
        else
            L1B.ACDC.SSH=L1BS.alt_sat-L1B.ACDC.tracker_range;
        end
        
        
        %---update the buffer: might be used for computation of initial params in a
        %sliding window approach
        geo_param_buffer(i_surf_stacked).epoch=L1B.ACDC.epoch;
        geo_param_buffer(i_surf_stacked).SSH=L1B.ACDC.SSH;
        if cnf_p_ACDC.rou_flag
            geo_param_buffer(i_surf_stacked).rou=estimates(2);
        else
            geo_param_buffer(i_surf_stacked).Hs=L1B.ACDC.Hs;
        end
        geo_param_buffer(i_surf_stacked).Pu=L1B.ACDC.Pu;
        geo_param_buffer(i_surf_stacked).error_epoch_ACDC=estimates(6);
    end
    %for the iterative approach
    Hs=L1B.ACDC.Hs;
    epoch=L1B.ACDC.epoch;    
%     disp(L1B.ACDC.Hs)
%     disp(L1B.ACDC.SSH)
end




%% --------------- UPDATING THE INITIAL_PARAMS ----------------------------
%--------------------------------------------------------------------------
if cnf_p_ACDC.initial_param_fit_feedback_flag
    % Amplitude fitting
    fit_params_ini(3) = 1;
    
    % initial epoch fitting
    fit_params_ini(1)  = estimates(1);
    
    % initial differential epoch
   
    
    
    % initial SWH/roughness seed    
    if i_surf_stacked~=1
        if cnf_p_ACDC.rou_flag
            average_total=nanmean([geo_param_buffer(1:i_surf_stacked-1).rou]);
            std_total=nanstd([geo_param_buffer(1:i_surf_stacked-1).rou]);
        else
            average_total=nanmean([geo_param_buffer(1:i_surf_stacked-1).Hs]./4);
            std_total=nanstd([geo_param_buffer(1:i_surf_stacked-1).Hs]./4);
        end
        
        average_total_error_epoch_ACDC=nanmean([geo_param_buffer(1:i_surf_stacked-1).error_epoch_ACDC]);
        std_total_error_epoch_ACDC=nanstd([geo_param_buffer(1:i_surf_stacked-1).error_epoch_ACDC]);
        
        average_total_SSH_ACDC=nanmean([geo_param_buffer(1:i_surf_stacked-1).SSH]);
        std_total_SSH_ACDC=nanstd([geo_param_buffer(1:i_surf_stacked-1).SSH]);
        
        if cnf_p_ACDC.ini_Hs_rou_sliding_win_opt==1
            %use sliding window to estimate the next initial seed
            init_sample=max(i_surf_stacked-(cnf_p_ACDC.ini_Hs_rou_sliding_win_size-1),1);
            if cnf_p_ACDC.rou_flag
                average_sliding=nanmean([geo_param_buffer(init_sample:i_surf_stacked-1).rou]);
                std_sliding=nanstd([geo_param_buffer(init_sample:i_surf_stacked-1).rou]);
            else
                average_sliding=nanmean([geo_param_buffer(init_sample:i_surf_stacked-1).Hs]./4);
                std_sliding=nanstd([geo_param_buffer(init_sample:i_surf_stacked-1).Hs]./4);
            end
                        
            %to discard effects of having a set of different
            %consecutive non-valid SWH/rou estimates
            if (average_sliding < average_total-cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total) || (average_sliding > average_total+cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total)
                %take the average from the very beginning
                fit_params_ini(2)=average_total;
            else
                current_value=L1B.ACDC.Hs/4;
                if ((current_value >= average_sliding-cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding) && (current_value <= average_sliding+cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding))
                    %if within the limits +- std of sliding include
                    %them in average
                    if cnf_p_ACDC.rou_flag
                        fit_params_ini(2)=nanmean([geo_param_buffer(init_sample:i_surf_stacked).rou]);
                    else
                        fit_params_ini(2)=nanmean([geo_param_buffer(init_sample:i_surf_stacked).Hs])/4.0;
                    end
                else
                    %if not wihtin limits +- std sliding take the
                    %sliding average
                    fit_params_ini(2)=average_sliding;
                end
            end
            
            %----------- Error epoch ACDC ---------------------------------
            average_sliding_error_epoch_ACDC=nanmean([geo_param_buffer(init_sample:i_surf_stacked-1).error_epoch_ACDC]);
            std_sliding_error_epoch_ACDC=nanstd([geo_param_buffer(init_sample:i_surf_stacked-1).error_epoch_ACDC]);
            
            if (average_sliding_error_epoch_ACDC < average_total_error_epoch_ACDC-cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total_error_epoch_ACDC) ||...
                    (average_sliding_error_epoch_ACDC > average_total_error_epoch_ACDC+cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total_error_epoch_ACDC)
                %take the average from the very beginning
                fit_params_ini(4)=average_total_error_epoch_ACDC;
            else
                current_value=estimates(6);
                if ((current_value >= average_sliding_error_epoch_ACDC-cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding_error_epoch_ACDC) && ...
                        (current_value <= average_sliding_error_epoch_ACDC+cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding_error_epoch_ACDC))
                    %if within the limits +- std of sliding include
                    %them in average
                    fit_params_ini(4)=nanmean([geo_param_buffer(init_sample:i_surf_stacked).error_epoch_ACDC]);
                else
                    %if not wihtin limits +- std sliding take the
                    %sliding average
                    fit_params_ini(4)=average_sliding_error_epoch_ACDC;
                end
            end
            
            %----------- SSH ---------------------------------
            average_sliding_SSH_ACDC=nanmean([geo_param_buffer(init_sample:i_surf_stacked-1).SSH]);
            std_sliding_SSH_ACDC=nanstd([geo_param_buffer(init_sample:i_surf_stacked-1).SSH]);
            
            if (average_sliding_SSH_ACDC < average_total_SSH_ACDC-cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total_SSH_ACDC) ||...
                    (average_sliding_SSH_ACDC > average_total_SSH_ACDC+cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total_SSH_ACDC)
                %take the average from the very beginning
                fit_params_ini(5)=average_total_SSH_ACDC;
            else
                current_value=L1B.ACDC.SSH;
                if ((current_value >= average_sliding_SSH_ACDC-cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding_SSH_ACDC) && ...
                        (current_value <= average_sliding_SSH_ACDC+cnf_p_ACDC.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding_SSH_ACDC))
                    %if within the limits +- std of sliding include
                    %them in average
                    fit_params_ini(5)=nanmean([geo_param_buffer(init_sample:i_surf_stacked).SSH]);
                else
                    %if not wihtin limits +- std sliding take the
                    %sliding average
                    fit_params_ini(5)=average_sliding_SSH_ACDC;
                end
            end
            
            
            
        else
            accumulated=accumulated+estimates(2)/4; %sigmas
            fit_params_ini(2) = (accumulated/i_surf_fitted); %sigmas
            accumulated_error_epoch=accumulated_error_epoch+estimates(6);
            fit_params_ini(4) = (accumulated_error_epoch/i_surf_fitted);
            accumulated_SSH=accumulated_SSH+L1B.ACDC.SSH;
            fit_params_ini(5) = (accumulated_SSH/i_surf_fitted);
        end
    else
        fit_params_ini(2)=estimates(2)/4; %sigmas
        fit_params_ini(4) = estimates(6);
        fit_params_ini(5) = L1B.ACDC.SSH;
    end
end


%------------------- verbose ----------------------------------------------
if cnf_p_ACDC.plot_fits_flag && ((mod(i_surf_stacked,cnf_p_ACDC.plot_fits_downsampling)==0) || i_surf_stacked==1)
    resultPath=files.resultPath;
    switch mission
        case 'CR2'
            filename = files.sph.product_info.product_id(20:20+30);
            if optional_ext_file_flag
                filename=strcat(filename,file_ext_string);
            end
    end
%     max_image=max([max(10*log10(stack_before_SRC(:))),max(10*log10(stack(:))),max(10*log10(stack_ac_before_smooth(:))),max(10*log10(stack_ac(:)))]);
%     f1=figure('NumberTitle','off','Name','Comparison Stack before and after ACDC');
%     title(strcat('Comparison Stack before and after ACDC #',{' '},num2str(i_surf_fitted-1)))
%     % subplot(2,2,1);
%     % surfc(k,l,10*log10(stack_before_SRC),'LineStyle','none','FaceColor','flat');
%     % title('Stack before SRC'); xlabel('Range index'); ylabel('Doppler index');
%     % colormap('jet'); c=colorbar; ylabel(c,'[dB]'); caxis([max_image-40.0, max_image]); view(2); grid off;
%     
%     subplot(2,2,1);
%     surfc(k,l,10*log10(stack),'LineStyle','none','FaceColor','flat');
%     title('Stack after SRC'); xlabel('Range index'); ylabel('Doppler index');
%     colormap('jet'); c=colorbar; ylabel(c,'[dB]'); caxis([max_image-40.0, max_image]); view(2); grid off;
%     
%     subplot(2,2,2);
%     surfc(k,l,10*log10(stack_ac_before_smooth),'LineStyle','none','FaceColor','flat');
%     title('Stack after AC'); xlabel('Range index'); ylabel('Doppler index');
%     colormap('jet'); c=colorbar; ylabel(c,'[dB]'); caxis([max_image-40.0, max_image]); view(2); grid off;
%     
%     subplot(2,2,3);
%     surfc(k,l,10*log10(stack_ac_smooth),'LineStyle','none','FaceColor','flat');
%     title('Stack after AC + Smoothing'); xlabel('Range index'); ylabel('Doppler index');
%     colormap('jet'); c=colorbar; ylabel(c,'[dB]'); caxis([max_image-40.0, max_image]); view(2); grid off;
%     
%     subplot(2,2,4);
%     surfc(k_scaled,l,10*log10(stack_ac_smooth),'LineStyle','none','FaceColor','flat');
%     title('Stack after ACDC'); xlabel('Range index'); ylabel('Doppler index');
%     colormap('jet'); c=colorbar; ylabel(c,'[dB]'); caxis([max_image-40.0, max_image]); view(2); grid off;
%     file_name=[resultPath,'plots/',filename,'_stack_comparison_',num2str(i_surf_stacked)];
% %    savefig(f1,[file_name,'.fig'])%
%     print('-dpng ',[file_name,'.png']);    
%     close(f1)
   
    f1=figure('NumberTitle','off','Name','Conventional Retracker fitting');        
    plot(L1B.ACDC.ml_wav_preliminary/max(L1B.ACDC.ml_wav_preliminary),'-b'); 
    hold on; plot(L1B.ACDC.ml_wav_fitted_preliminary,'-r') 
    xlabel('Range bins'); ylabel('Norm. amplitude'); 
    legend('Conventional Waveform','Conventional fitting'); 
    title(['Conventional retracker over surface #',{' '},num2str(i_surf_fitted-1)])
    annotation('textbox', [0.725 0.25 0.18 0.18],...
        'String',{['SSH= ', num2str(L1B.ACDC.SSH_preliminary,4), ' [m]'],...
        ['SWH = ' ,num2str(L1B.ACDC.Hs_preliminary,4), ' [m]']},...
        'FontSize',18,...
        'FontName','Helvetica',...
        'BackgroundColor',[1 1 1]);
    file_name=[resultPath,'plots/',filename,'_conventional_fitting_',num2str(i_surf_stacked)];
    %savefig(f1,[file_name,'.fig'])%
    print('-dpng ',[file_name,'.png']);
    close(f1);
    
    f1=figure('NumberTitle','off','Name','ACDC Retracker fitting');        
    plot(L1B.ACDC.waveform/max(L1B.ACDC.waveform),'-b'); 
    hold on; plot(L1B.ACDC.ml_wav_fitted_ACDC,'-r') 
    xlabel('Range bins'); ylabel('Norm. amplitude'); 
    legend('ACDC waveform','ACDC fitting'); 
    title(['ACDC retracker over surface #',{' '},num2str(i_surf_fitted-1)])
        annotation('textbox', [0.725 0.25 0.18 0.18],...
        'String',{['SSH= ', num2str(L1B.ACDC.SSH,4), ' [m]'],...
        ['SWH = ' ,num2str(L1B.ACDC.Hs,4), ' [m]']},...
        'FontSize',18,...
        'FontName','Helvetica',...
        'BackgroundColor',[1 1 1]);
    file_name=[resultPath,'plots/',filename,'_ACDC_fitting_',num2str(i_surf_stacked)];
    %savefig(f1,[file_name,'.fig'])%
    print('-dpng ',[file_name,'.png']);
    close(f1);
    
end

end %end of function
    
% switch ACDC_ref_beam_method_est_cnf
%     case 'max'
%         %Based on the averaged value across-track get the beam with maximum
%         %averaged values
% 
%     case 'nadir'
%         % Based on the look angle information: look angle close to zero
%         [~,ref_beam]=min(abs(L1BS.look_ang_surf));
% end

