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
% FITFUNCTIONS: Function that computes the deltas in the x-axis (deltaT)
% and the y-axis (deltaR) that minimize the RMS between the two input
% signals (rtheo and rdata)
% 
% Calling
%   [deltaR, deltaT] = fitfunctions(rtheo , rdata)
%
% Inputs
%   rtheo : theoretical range
%   rdata : measured range
%
% Outputs
%   deltaR: delta range that minimizes the rms between the measured and the theoretical ranges
%   deltaT: delta time that minimizes the rms between the measured and the theoretical datation
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
% Reviewer: Monica Roca / isardSAT
%
%
%
% $Id: fitfunctions.m 139 2005-12-20 07:48:18Z merche $
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [deltaR, deltaT] = fitfunctions(rtheo , rdata)

  [coef]=  fminsearch_stackcal('fitfunc', [0 0] , [], rtheo, rdata)
  deltaT = coef(1);
  deltaR = coef(2);
cd ..\stack_data

