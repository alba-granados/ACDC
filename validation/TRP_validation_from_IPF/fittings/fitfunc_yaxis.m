% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
% 
% This code implements the algorithm as described in the
% ISARD_ESA_CR2_TRP_CAL_DPM_030 2.b of 26/05/2011
%
% ---------------------------------------------------------
% FITFUNC_YAXIS: function that computes the RMS with the delta_AOA (x0) coefficient from the input
% 
% Calling
%   E = fitfunc_yaxis (x0,AOAtheo,AOAdata)
%
% Inputs
%   x0 : deltaAOA
%   AOAtheo : theoretical angle
%   AOAdata : measured angle
%
% Output
%   E : rms
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
% Reviewer: Monica Roca / isardSAT
%
%
%
% $Id: fitfunc_yaxis.m 139 2005-12-20 07:48:18Z merche $
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = fitfunc_yaxis(x0, AOAtheo,AOAdata)
    sum = 0;
    for j=1:length(AOAtheo)  
        sum=sum+((AOAdata(j)+x0)- AOAtheo(j)).^2;    
    end
    E=sum;


