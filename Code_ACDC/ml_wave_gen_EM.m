function [ml_waveform,stack,stack_real] = ml_wave_gen_EM(x,fit_params,nf_p,cnf_p,looks,func_f0,func_f1)%varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L.
% -------------------------------------------------------------------------
% Generation of multi-looked waveform in IEEE TGRS "SAR Altimeter Backscattered Waveform Model"
% DOI:10.1109/TGRS.2014.23330423
% -------------------------------------------------------------------------
%
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Cristina Martin-Puig / isardSAT
%
% Last revision:    Cristina Martin-Puig / isardSAT V9 4/7/2014
%
% This software is built with internal funding
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       x           =   Range bin index
%       fit_p       =   Parameters to be fitted
%       nf_p        =   parameters of need for the waveform generation, but not to be
%                       fitted
%       cnf_p       =   configuration parameters of L2 processor
% OUTPUT:
%       ml_waveform =   multi-looked waveform
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - sl_wave_gen_EM: in charge of genrating the single look waveform based on
% the model proposed by Chris in IEEE TGRS "SAR Altimeter Backscattered Waveform Model"
% DOI:10.1109/TGRS.2014.23330423
% -
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Versions control:
% Based on the V9 version defined by Cristina Martin-Puig and updated
% according to the findings during the Sentinel-6 analysis (look indexation)


if nargin < 7
    func_f1=0; % initialize altitude to 0;
end
if nargin <6
    func_f0=0; % initialize altitude to 0;
end

%--------------------------------------------------------------------------
%--------------------- Stack modeling -------------------------------------
%--------------------------------------------------------------------------
n_looks=length(looks);
n_samples=length(x);
stack       =   zeros(n_looks,n_samples);
% for l = 1:n_looks
%     %nf_p.ThN = 0;
%     %modified by EM 01.03.2016 / modified by EM 15.03.2016
%     if cnf_p.lut_flag
%         switch cnf_p.power_wfm_model
%             case 'simple'
%                 stack(l,:)   = sl_wave_gen_EM(x,looks(l),fit_params,nf_p,cnf_p,func_f0);
%             case 'complete'
%                 stack(l,:)   = sl_wave_gen_EM(x,looks(l),fit_params,nf_p,cnf_p,func_f0,func_f1);
%         end
%     else
%         stack(l,:)   = sl_wave_gen_EM(x,looks(l),fit_params,nf_p,cnf_p);
%     end        
%     
% end
x_matrix=ones(n_looks,1)*x;
looks_matrix=looks*ones(1,n_samples);
if cnf_p.lut_flag
    switch cnf_p.power_wfm_model
        case 'simple'
            stack   = stack_gen(x_matrix,looks_matrix,fit_params,nf_p,cnf_p,func_f0);
        case 'complete'
            stack   = stack_gen(x_matrix,looks_matrix,fit_params,nf_p,cnf_p,func_f0,func_f1);
    end
else
    stack   = stack_gen(x_matrix,looks_matrix,fit_params,nf_p,cnf_p);
end
    

switch cnf_p.multilook_option
    case 'Cris_old'
        stack= fit_params(3).*stack+(nf_p.ThN).'*ones(1,n_samples);
end

% -------------------------------------------------------------------------
% ----------------------- Stack masking -----------------------------------
% -------------------------------------------------------------------------
stack_real     =   stack.*nf_p.Doppler_mask;
ml_waveform = nanmean(stack_real(1:length(looks),:),1);

%patch EM: 14.09.2016
%due to the way the stack mask vector is constructed and saved in the L1B
%product the last sample is always forced to zero (?¿) and so this implies
%that as it is saved in an integer when reading samples 255 and 256 will be
%set to zero differently from what is performed at level 1B
idx_nans=isnan(ml_waveform);
if any(idx_nans)
    idx_ref=find(~idx_nans,1,'last');
    %force the NaN samples to the same value as per last non-NaN value
    ml_waveform(idx_nans)=ml_waveform(idx_ref);
end


% -------------------------------------------------------------------------
% ----------------------- Multilooking ------------------------------------
% -------------------------------------------------------------------------
switch cnf_p.multilook_option
    case 'NormML_plus_noise'
        max_ML_signal=max(ml_waveform);        
        ml_waveform = (ml_waveform/max_ML_signal+nf_p.ThN(1));
        ml_waveform = ml_waveform/max(ml_waveform);

end


end

