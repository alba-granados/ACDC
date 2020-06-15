%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop
% This code implements the dilation compensation step / re-scaling of the ACDC stack
% generation based on paper of Chris et al. 2015
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
% Last revision:    Eduard Makhoul / isardSAT V1 29/09/2016
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

function [k_scaled]=dilation_compensation(k,g,nf_p,SWH,look_indexation_ref)

%% ------------Computation of g0 ------------------------------------------
sigmaz  =   SWH/4.0; 
l_0=look_indexation_ref;
g_0       =   sqrt(2*nf_p.alphag_a*nf_p.alphag_r./(nf_p.alphag_a + nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * l_0.^2 + 2 * nf_p.alphag_a * nf_p.alphag_r * (sigmaz/nf_p.Lz)^2 ));        

%% ---------- Scaling of the range bin index ------------------------------
k_scaled=g.*k/g_0;



end


