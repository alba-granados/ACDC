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
% FITFUNCTIONS_YAXIS: Function that computes the deltas in the y-axis that minimize the RMS between the two input
% signals (AOAtheo and AOAdata)
% 
% Calling
%   delta_AOA = fitfunctions_yaxis(AOAtheo , AOAdata)
%
% Inputs
%   AOAtheo : theoretical angle
%   AOAdata : measured angle
%
% Outputs
%   delta_AOA: delta angle that minimizes the rms between the measured and the theoretical angles
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
% Reviewer: Monica Roca / isardSAT
%
% Last revision: Josep Montolio / Pildo Labs (19/10/09)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta_Y = fitfunctions_yaxis(AOAtheo , AOAdata)

  delta_Y=  fminsearch_stackcal('fitfunc_yaxis', 0 , [], AOAtheo, AOAdata);



