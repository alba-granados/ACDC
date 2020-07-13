%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop
% This code implements the fitting to a single look waveform model with
% first order approximation after ACDC A*f0()
%
% ---------------------------------------------------------
% Objective: Perform the amplitude compensation of the stack
% 
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Roger Escol� / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (02/06/2015)
% 
% 
% Last revision:    Eduard Makhoul / isardSAT V1 30/09/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -waveform = input waveform to be fitted
%       -k = range vector indexation - estimated epoch bin in a matricial
%       notation
%       -fit_params_ini = initial parameters [error_kappa,Hs,Pu]
%       -nf_p = non-fitting parameters strucrture
%       -SWH = estimated SWH from the preliminary fitting
%       
% OUTPUT:
%       -estimates = vector of parameter estimates [error k estimation, Hs (m), Pu (dB)]
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

function [estimates,ml_wav]=fit_sl_f0_model_ACDC(waveform,k,nf_p,cnf_p_ACDC,fit_params_ini,look_index_ref,i_surf_stacked,func_f0)
if nargin < 8
    func_f0=0.0;
end


%% -------------------- NORMALIZATION OF DATA -----------------------------
%--------------------------------------------------------------------------
max_waveform=max(waveform);
fit_data=waveform./max_waveform;

%fit_params_ini = initial parameters [epoch,Hs,Pu]

%-------------------- PROVIDE initial fittings conversion -----------------
sigmaz = fit_params_ini(2)/4;
fit_params_ini(2) = sqrt(2*nf_p.alphag_a*nf_p.alphag_r./(nf_p.alphag_a + nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * look_index_ref.^2 + 2 * nf_p.alphag_a * nf_p.alphag_r * (sigmaz/nf_p.Lz)^2 ));  %define the gl parameter to use in fitting

%% -------------------- PRE-PROCESSING & NOISE ESTIMATION -----------------
% -------------------------------------------------------------------------
% Pre-processing
% -------------------------------------------------------------------------
if cnf_p_ACDC.pre_processing
    %used to estimate the initial value of the epoch 
    %& eventually if desired the noise floor
    [peak_pow,idx_max_peak]=max(fit_data);
%     plot(fit_data);
    idx_leading_edge=find(find(fit_data<=cnf_p_ACDC.percent_leading_edge/100.0*peak_pow)<idx_max_peak, 1, 'last' );
    if isempty(idx_leading_edge)
        %if there is no leading edge or the waveform has displaced that
        %much from the window to the left select the peak as leading
        %edge
        idx_leading_edge=idx_max_peak;
    end
    %fit_params_ini(1)=k(idx_leading_edge);
    
    %noise floor estimation
    if cnf_p_ACDC.fit_noise
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
                    ThN        =   mean(fit_data(idx_noise));
                    clear idx_noise tempo_noise_thr;
                else
                    %take an average mean value of the previous Thermal
                    %noise estimation and keep it as reference for this
                    %record
                    if i_surf_stacked~=1
                        ThN         = mean(nf_p.ThN);
                    else
                        %!!!need to be reviewed
                        %peak at the beginning of the window and so no
                        %space left for estimation of the noise floor
                        ThN         = min(fit_data);
                    end
                end
        end
    end
end
%--------------------- Using the same noise window for all surfaces -------
if cnf_p_ACDC.Thn_flag && ~cnf_p_ACDC.fit_noise
    ThN        =   mean(fit_data(cnf_p_ACDC.Thn_w_first:cnf_p_ACDC.Thn_w_first+cnf_p_ACDC.Thn_w_width-1));
end


%% ------------------------ FITTING ---------------------------------------
% -------------------------------------------------------------------------
% DEFINE FITTING MODEL & CALL FITTING ROUTINE 
% -------------------------------------------------------------------------                        
% -------------------------------------------------------------
% --------- DEFINE FITTING FUNCTION MODEL ---------------------
% -------------------------------------------------------------
%[waveform_ml_ACDC]=sl_f0_model_ACDC(k,param,ThN,cnf_p_ACDC,func_f0)
if cnf_p_ACDC.lut_flag
    switch cnf_p_ACDC.fitting_fun_type
        case 'lsq'
            mpfun               =   @(fit_params,x)sl_f0_model_ACDC(x, fit_params, ThN, cnf_p_ACDC,func_f0);
        case 'fmin'
            fminfun             =   @(fit_params)sum((sl_f0_model_ACDC(k, fit_params, ThN, cnf_p_ACDC,func_f0)-fit_data).^2);
    end
    
else
    switch cnf_p_ACDC.fitting_fun_type
        case 'lsq'
            mpfun               =   @(fit_params,x)sl_f0_model_ACDC(x, fit_params, ThN, cnf_p_ACDC);
        case 'fmin'
            fminfun             =   @(fit_params)sum((sl_f0_model_ACDC(k, fit_params, ThN, cnf_p_ACDC)-fit_data).^2);
    end
end

% -------------------------------------------------------------------------
% --------- RUN THE MINIMIZATION/FITTING PROCEDURE ------------------------
% -------------------------------------------------------------------------
% disp('input of lsqcurvefit:');
% disp(fit_params_ini)
% figure(); plot(fit_data);

bound = sqrt(2 * nf_p.alphag_a * nf_p.alphag_r/(nf_p.alphag_a + nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * look_index_ref.^2));
cnf_p_ACDC.fitting_options_ACDC_lb(2) = -bound;
cnf_p_ACDC.fitting_options_ACDC_ub(2) = bound;
switch cnf_p_ACDC.fitting_fun_type
    case 'lsq'
        [estimates,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini,k,fit_data,...
                                                   cnf_p_ACDC.fitting_options_ACDC_lb,cnf_p_ACDC.fitting_options_ACDC_ub,cnf_p_ACDC.fitting_options);
    case 'fmin'
        [estimates,~,flag]     =   fminsearchbnd (fminfun,fit_params_ini,zeros(1,length(fit_params_ini)),...
                                                   cnf_p_ACDC.fitting_options_ACDC_lb,cnf_p_ACDC.fitting_options);
end
% disp('output of lsqcurvefit:')
% disp(estimates)
% disp('() > 0?:');
% disp(nf_p.Lz^2/(2 * nf_p.alphag_a * nf_p.alphag_r)*(2*nf_p.alphag_a*nf_p.alphag_r./(estimates(2)^2)-nf_p.alphag_a-nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * look_index_ref.^2));
%% --------- DEFINE THE OUTPUT PARAMETERS ---------------------------------
%--------------------------------------------------------------------------
%---------- computation of the correlation coefficient --------------------
if cnf_p_ACDC.lut_flag
    ml_wav=sl_f0_model_ACDC(k, estimates, ThN, cnf_p_ACDC,func_f0);
else
    ml_wav=sl_f0_model_ACDC(k, estimates, ThN, cnf_p_ACDC);
end
correlation_fit         =   corrcoef(fit_data,ml_wav);
% estimates(1): refers to the error in epoch estimation
% estimates(2): significant waveheight to be extracted from gl estimated
estimates(2)=4*sqrt(nf_p.Lz^2/(2 * nf_p.alphag_a * nf_p.alphag_r)*(2*nf_p.alphag_a*nf_p.alphag_r./(estimates(2)^2)-nf_p.alphag_a-nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * look_index_ref.^2));
% estimates(3): Pu amplitude
estimates(3)=10*log10(estimates(3)*max_waveform); %dB
estimates(4) = correlation_fit(1,2)*100;

%fitting flag
estimates(5)=flag;

end

