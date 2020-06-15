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
% FITFUNC_XAXIS: function that computes the RMS with the delta_AOA (x0) coefficient from the input
% 
% Calling
%   E = fitfunc_xaxis (x0,AOAtheo,AOAdata)
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
% $Id: fitfunc_xaxis.m 139 2005-12-20 07:48:18Z merche $
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function E = fitfunc_xaxis(x0, theo, meas)
    sum = 0;
    x = abs (round(x0(1)));
    [~,minPos]=min(theo);
    X_margin=1000;
    for j=(minPos-X_margin):(minPos+X_margin)  
       
%         if (j + x <= length(meas)) 
            sum=sum+((meas(j+x))- theo(j)).^2;    
%         end
    end
%     figure;
%     hold on; p2=plot((minPos-X_margin):(minPos+X_margin),theo((minPos-X_margin):(minPos+X_margin)),'g');
%     set(p2,{'DisplayName'},{sprintf('theorethical')})
%     hold on; p=plot((minPos-X_margin):(minPos+X_margin),meas((minPos-X_margin+x):(minPos+X_margin+x)));
%     set(p,{'DisplayName'},{sprintf('measured + %d = %d',x0,sum)})

    E=sum;


