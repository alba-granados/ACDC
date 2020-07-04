%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop
% This code implements the fitting of simplified Chris Ray model f0 to
% waveforms in the stack (with no multilook) to obtain an average estimation of the 
% Hs and koff
%
% ---------------------------------------------------------
% Objective: Provide a flattened version of the stack a la ACDC and the
% corresponding multilooked ACDC waveform
% 
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Roger Escol� / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (02/06/2015)
% 
% 
% Last revision:    Eduard Makhoul / isardSAT V1 29/09/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -stack = input masked stack after geometry corrections and range
%       -stack_mask = input stack mask 
%       -looks = looks indexation 
%       -nf_p = non-fitting parameters strucrture
%       -cnf_p = configuration parameters for fitting
%       
% OUTPUT:
%       -SWH        =   preliminary estimate of SWH
%       -epoch      =   preliminary estimate of the epoch
%       -amplitude  =   fitted amplitude
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - Depending on the stack being analyzed (first surfaces) is not always
% possible to have all the left and right beams: truncate to the available
% possible ones avoiding those set
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0:
% V1.1: 

function [Hs, epoch, Pu,corr_coeff, fit_params_ini, nf_p, flag_exit,ml_wav,ml_wav_fitted]  = stack_fitting(stack,stack_mask,looks,range_index,...
                                                                          nf_p,cnf_p_ACDC,fit_params_ini,i_surf_stacked,...
                                                                          func_f0,func_f1)
if nargin < 10
    func_f1=0; % initialize altitude to 0;
end
if nargin < 9
    func_f0=0; % initialize altitude to 0;
end

global avoid_beams_mask_allzeros_cnf
flag_exit=0;
%% ------------------- Formation of the Multilooked waveform --------------
%--------------------------------------------------------------------------
%Since we don't have an average over several bursts we may consider the
%multilooked waveform
% If we perform the fitting over only one waveform is highly affected by
% speckle
if avoid_beams_mask_allzeros_cnf % alba: is this already computed in multilooking.m, see L1B.start_beam, .stop_beam?
    %mask: beams where all the samples are set to 0 are non-contributing to the
    %ML
    idx_rows_int=any(stack_mask,2);
else
    idx_rows_int=1:1:nf_p.Neff;
end
fit_data=nanmean(stack(idx_rows_int,:),1);
ml_wav=fit_data;
%for the fitting
max_fit_data=nanmax(fit_data);
fit_data=fit_data./max_fit_data;

%check whether all samples waveform are NaN due to the masking and or due to fact
%that noisy beams removal filters all them
if ~any(~isnan(fit_data))
    flag_exit=-1;
    Hs=NaN;
    epoch=NaN;
    a=NaN;
    return;
end

%-------------------------- Load the mask into parameter ------------------
nf_p.Doppler_mask=stack_mask;


%% -------------------- PRE-PROCESSING & NOISE ESTIMATION -----------------
% -------------------------------------------------------------------------
% Pre-processing
% -------------------------------------------------------------------------
if cnf_p_ACDC.pre_processing
    %used to estimate the initial value of the epoch 
    %& eventually if desired the noise floor
    fprintf('estimating initial value of the epoch...\n')
    [peak_pow,idx_max_peak]=max(fit_data);
    idx_leading_edge=find(find(fit_data<=cnf_p_ACDC.percent_leading_edge/100.0*peak_pow)<idx_max_peak, 1, 'last' );
    if isempty(idx_leading_edge)
        %if there is no leading edge or the waveform has displaced that
        %much from the window to the left select the peak as leading
        %edge
        idx_leading_edge=idx_max_peak;
    end
    fit_params_ini(1)=range_index(idx_leading_edge);
    
    %noise floor estimation
    if cnf_p_ACDC.fit_noise
        fprintf('estimating thermal noise...\n')
        switch cnf_p_ACDC.Thn_method
            case 'ML'
                idx_noise=[];
                temp_noise_thr=cnf_p_ACDC.threshold_noise;
                iter_noise=1;
                while isempty(idx_noise) && iter_noise<cnf_p_ACDC.max_iter_noise
                    idx_noise=find(abs(diff(fit_data))<=temp_noise_thr);
                    idx_noise=idx_noise(idx_noise<idx_leading_edge);
                    temp_noise_thr=temp_noise_thr*1.5;
                    iter_noise=iter_noise+1;
                end
                if iter_noise<cnf_p_ACDC.max_iter_noise
                    nf_p.ThN        =   mean(fit_data(idx_noise))*ones(1,length(looks));
                    clear idx_noise tempo_noise_thr;
                else
                    %take an average mean value of the previous Thermal
                    %noise estimation and keep it as reference for this
                    %record
                    if i_surf_stacked~=1
                        nf_p.ThN         = mean(nf_p.ThN)*ones(1,length(looks));
                    else
                        %!!!need to be reviewed
                        %peak at the beginning of the window and so no
                        %space left for estimation of the noise floor
                        nf_p.ThN         = min(fit_data)*ones(1,length(looks));
                    end
                end
        end
    end
end
%--------------------- Using the same noise window for all surfaces -------
if cnf_p_ACDC.Thn_flag && ~cnf_p_ACDC.fit_noise
    % added by EM 02/12/2015 generate a vector of thermal noise for the
    % different looks
    switch cnf_p_ACDC.Thn_method
        case 'ML'
            nf_p.ThN        =   mean(fit_data(cnf_p_ACDC.Thn_w_first:cnf_p_ACDC.Thn_w_first+cnf_p_ACDC.Thn_w_width-1))*ones(1,nf_p.Neff);
        case 'SL'
            dumm=squeeze(stack(1:nf_p.Neff,:));
            max_mean_dumm=max(nanmean(dumm,1));
            for i_beam=1:nf_p.Neff(m)
                nf_p.ThN(i_beam)        =   mean(dumm(i_beam,cnf_p_ACDC.Thn_w_first:cnf_p_ACDC.Thn_w_first+cnf_p_ACDC.Thn_w_width-1)/max_mean_dumm,'omitnan');
            end
            
            clear dumm;
    end
end


%% ------------------------ FITTING ---------------------------------------
% -------------------------------------------------------------------------
% DEFINE FITTING MODEL & CALL FITTING ROUTINE 
% -------------------------------------------------------------------------                        
% -------------------------------------------------------------
% --------- DEFINE FITTING FUNCTION MODEL ---------------------
% -------------------------------------------------------------
if cnf_p_ACDC.lut_flag
    fprintf('defining fitting fuctions model...\n')
    switch cnf_p_ACDC.power_wfm_model
        case 'simple'
            switch cnf_p_ACDC.fitting_fun_type
                case 'lsq'
                    mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_p, cnf_p_ACDC,looks,func_f0);
                case 'fmin'
                    fminfun             =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_p, cnf_p_ACDC,looks,func_f0)-fit_data).^2);
            end
            
        case 'complete'
            switch cnf_p_ACDC.fitting_fun_type
                case 'lsq'
                    mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_p, cnf_p_ACDC,looks,func_f0,func_f1);
                case 'fmin'
                    fminfun             =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_p, cnf_p_ACDC,looks,func_f0,func_f1)-fit_data).^2);
            end
    end
else
    switch cnf_p_ACDC.fitting_fun_type
        case 'lsq'
            mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_p, cnf_p_ACDC,looks);
        case 'fmin'
            fminfun             =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_p, cnf_p_ACDC,looks)-fit_data).^2);
    end
end

% -------------------------------------------------------------------------
% --------- RUN THE MINIMIZATION/FITTING PROCEDURE ------------------------
% -------------------------------------------------------------------------
switch cnf_p_ACDC.fitting_fun_type
    case 'lsq'
        fprintf('running minimization %s procedure...\n', cnf_p_ACDC.fitting_fun_type);
        [fit_params,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini,range_index,fit_data,...
                                                   cnf_p_ACDC.fitting_options_lb,cnf_p_ACDC.fitting_options_ub,cnf_p_ACDC.fitting_options);
    case 'fmin'
        [fit_params,~,flag]     =   fminsearchbnd (fminfun,fit_params_ini,zeros(1,length(fit_params_ini)),...
                                                   cnf_p_ACDC.fitting_options_ub,cnf_p_ACDC.fitting_options);
end


%% --------- DEFINE THE OUTPUT PARAMETERS ---------------------------------
%--------------------------------------------------------------------------
epoch    =   fit_params(1); % include again the IF mask contribution
%fit_params_ini(1)=epoch;
Hs      =   (abs(fit_params(2))*4); % sigma_z = Hs/4
if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
    Pu       =   10*log10(fit_params(3)*max_fit_data); %[dB]
else
    Pu       =   10*log10(max_fit_data);%[dB]
end

%--------- Pearson correlation coefficient --------------------------------
if cnf_p_ACDC.lut_flag
    switch cnf_p_ACDC.power_wfm_model
        case 'simple'
            ml_wav_fitted=ml_wave_gen_EM(range_index, fit_params, nf_p, cnf_p_ACDC,looks,func_f0);
            
        case 'complete'
            ml_wav_fitted=ml_wave_gen_EM(range_index, fit_params, nf_p, cnf_p_ACDC,looks,func_f0,func_f1);
    end
else
    ml_wav_fitted=ml_wave_gen_EM(range_index, fit_params, nf_p, cnf_p_ACDC,looks);
end
correlation_fit         =   corrcoef(fit_data,ml_wav_fitted);
corr_coeff              =   correlation_fit(1,2)*100;

end


