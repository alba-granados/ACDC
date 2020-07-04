%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop
% This code implements the generation of a single look waveform model with
% first order approximation after ACDC A*f0(g*k)
%
% ---------------------------------------------------------
% Objective:
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
%       -k     = kappa values (range indexes)
%       -param = [offset kappa, gl , Pu]
%       -cnf_p_ACDC = structure of processing parameters
%      OPTIONAL
%       -func_f0 = tabulated LUT of f0
%       
% OUTPUT:
%       -waveform_ml_ACDC = modeled waveform param(3)*f0(g(param(2))*(k-param(1))) 
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

function [waveform_ml_ACDC]=sl_f0_model_ACDC(k,param,ThN,cnf_p_ACDC,func_f0)
if nargin < 5
    func_f0=0;
end

gk=(k-param(1))*param(2);
waveform_ml_ACDC=0.*gk;
if cnf_p_ACDC.lut_flag    
    waveform_ml_ACDC(gk==0)=1.077900274770464; % 
    indexes_1=gk>=cnf_p_ACDC.LUT_ximin & gk<=cnf_p_ACDC.LUT_ximax;
    indexes_2=floor((gk(indexes_1)-cnf_p_ACDC.LUT_ximin)./cnf_p_ACDC.LUT_step)+1;
    indexes_3=gk>cnf_p_ACDC.LUT_ximax;
    waveform_ml_ACDC(indexes_1)=func_f0(indexes_2);
    waveform_ml_ACDC(indexes_3)=sqrt(pi./(2.0*gk(indexes_3))).*(1+3./(8*(gk(indexes_3)).^2)+105./(16*8*(gk(indexes_3)).^4));
else
    waveform_ml_ACDC(gk~=0)=pi/4.0*sqrt(abs(gk(gk~=0))).*(besseli(-1/4,1/4*(gk(gk~=0)).^2,1)+sign(gk(gk~=0)).*besseli(1/4,1/4*(gk(gk~=0)).^2,1)); 
    waveform_ml_ACDC(gk==0)=1.077900274770464;%2^(1/4)*gamma(5/4);    
end

waveform_ml_ACDC=param(3).*waveform_ml_ACDC+ThN;



end

