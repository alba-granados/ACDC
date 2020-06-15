%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop
% This code implements the retracking of the ACDC waveform after AC and DC
% have been applied: two step fitting
%
% ---------------------------------------------------------
% Objective: Perform the amplitude compensation of the stack
% 
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Roger Escolà / isardSAT
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
%       -l = looks indexation in a matricial notation
%       -k = range vector indexation - estimated epoch bin in a matricial
%       notation
%       -g = dilation term for the estimated SWH in matricial notation
%       -nf_p = non-fitting parameters strucrture
%       -SWH = estimated SWH from the preliminary fitting
%       
% OUTPUT:
%       -k_scaled = scaled range indexation (-epoch), including the 
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

function [waveform_ml_ACDC,k_ml_ACDC,estimates,flag_exit,ml_wav_fitted]=ACDC_retracking(stack_ac,...
                                                                k_scaled,nf_p,cnf_p_ACDC,fit_params_ini,...
                                                                look_index_ref,i_surf_stacked,func_f0)
if nargin < 8
    func_f0=0.0;
end

flag_exit=0;
global N_samples zp_fact_range_cnf

%% -------------------- FIRST LEVEL FITTING ACDC STACK --------------------
%--------------------------------------------------------------------------
% Re-order the AC stack
stack_1D=reshape(stack_ac,1,nf_p.Neff*(N_samples*zp_fact_range_cnf));
idx_finite=isfinite(stack_1D);
stack_1D=stack_1D(idx_finite);
N_samples_stack_1D=length(find(idx_finite));
% Re-order the scaled range
k_scaled_1D=reshape(k_scaled,1,nf_p.Neff*(N_samples*zp_fact_range_cnf));
k_scaled_1D=k_scaled_1D(idx_finite);
% Run the fitting
%[estimates,ml_wav_fitted]=fit_sl_f0_model_ACDC(stack_1D,k_scaled_1D,nf_p,cnf_p_ACDC,[0,fit_params_ini(2),1.0],look_index_ref,i_surf_stacked,func_f0);


%% -------------------- MULTILOOKING OF THE ACDC --------------------------
%--------------------------------------------------------------------------
kmin = floor(min(k_scaled_1D));
kmax = floor(max(k_scaled_1D));
switch lower(cnf_p_ACDC.weighting_win_type)
    case 'gaussian'
        k_ml_ACDC = (kmin:cnf_p_ACDC.weighting_win_width:kmax);
        N_samples_ACDC = length(k_ml_ACDC);
        matrix_weighting=exp(-(1.0*((ones(N_samples_ACDC,1)*k_scaled_1D)-((k_ml_ACDC.')*ones(1,N_samples_stack_1D)))).^2/(2*cnf_p_ACDC.weighting_win_width.^2));
        matrix_weighting=matrix_weighting./(matrix_weighting*ones(N_samples_stack_1D,1)*ones(1,N_samples_stack_1D));
        %waveform_ml_ACDC=((ones(N_samples_ACDC,1)*stack_1D)*matrix_weighting.').';
        waveform_ml_ACDC=(matrix_weighting*(stack_1D.')).';
    otherwise
        error('No valid weighting function')    
end


%% -------------------- FITTING OF THE MULTILOOKED ACDC -------------------
%--------------------------------------------------------------------------
[estimates,ml_wav_fitted]=fit_sl_f0_model_ACDC(waveform_ml_ACDC,k_ml_ACDC,nf_p,cnf_p_ACDC,[fit_params_ini(4),fit_params_ini(2),1.0],look_index_ref,i_surf_stacked,func_f0);

%differential error of the epoch
estimates(6)=estimates(1);

%epoch: add the estimated error over the initial epoch
estimates(1)=estimates(1)+fit_params_ini(1);



%----------- re-organize the waveforms (may not be original size)----------
dumm=waveform_ml_ACDC;
waveform_ml_ACDC=ones(1,N_samples*zp_fact_range_cnf)*NaN;
dumm2=ml_wav_fitted;
ml_wav_fitted=ones(1,N_samples*zp_fact_range_cnf)*NaN;
dumm3=k_ml_ACDC;
k_ml_ACDC=ones(1,N_samples*zp_fact_range_cnf)*NaN;
if N_samples_ACDC > N_samples*zp_fact_range_cnf
   waveform_ml_ACDC=dumm(1:N_samples*zp_fact_range_cnf);
   ml_wav_fitted=dumm2(1:N_samples*zp_fact_range_cnf);
   k_ml_ACDC=dumm3(1:N_samples*zp_fact_range_cnf);
else
   waveform_ml_ACDC(1:N_samples_ACDC)=dumm(1:N_samples_ACDC);
   ml_wav_fitted(1:N_samples_ACDC)=dumm2(1:N_samples_ACDC);
   k_ml_ACDC(1:N_samples_ACDC)=dumm3(1:N_samples_ACDC);
end


end