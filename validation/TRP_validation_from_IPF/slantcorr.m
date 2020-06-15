% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% -----------------------------------------------  ---------
%
% CryoSat 2 calibration over transponders
% 
% This code implements the algorithm as described in the
% ISARD_ESA_CR2_TRP_CAL_DPM_030 2.b of 26/05/2011
%
% ---------------------------------------------------------

% SLANTCORR: computes the slant correction
% 
% Calling
%   range_out = slantcorr(wvf_stack_ini, range_in)
%
% Inputs
%    wvf_stack_ini: stack corresponding to the TRP
%    range_in: range with the slant correction applied
%
% Outputs
%    range_out: range with the slant correction removed
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
%           Josep Montolio / Pildo Labs
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Mònica Roca / isardSAT (26/05/11)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function range_out = slantcorr(wvf_stack_ini, range_in, nbeams, start, final)
  

  R_dop = wvf_stack_ini.doppler_corr(start:nbeams-final)* 1000;
  R_dop = R_dop';
  R_dop =0;
%   disp('doppler corr set to 0 to datation issue');
  
  R_sla = wvf_stack_ini.slant_range(start:nbeams-final) * 1000;    % millimeters
  R_sla = R_sla';
  
  range_out = range_in + R_sla + R_dop;
