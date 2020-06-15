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
% FITFUNC: function that computes the RMS with the deltaT (x0(1)) and deltaR(x0(2))
% coefficients from the input
% 
% Calling
%   E = fitfunc(x0,rtheo,rdata3)
%
% Inputs
%   x0 : deltaR and deltaT 
%   rtheo : theoretical range
%   rdata3 : measured range
%
% Output
%   E : rms
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
% Reviewer: Monica Roca / isardSAT
%
% $Id: fitfunc.m 139 2005-12-20 07:48:18Z merche $
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function E = fitfunc(x0,rtheo,rdata3)
sum = 0;
x0   % x0(1) = deltaT 
     % x0(2) = deltaR 

 for j=1: length(rdata3)
            x = abs(round(x0(1)));
            if (j + x <= length(rdata3)) 
                sum=sum+(((rtheo(j+x))+ x0(2))- rdata3(j)).^2;    
            end    
            if (j + x > length(rdata3)), break, end
 end
 E=sum;


