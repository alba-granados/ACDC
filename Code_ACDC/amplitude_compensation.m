%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop
% This code implements the amplitude compensation step of the ACDC stack
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
%       -stack = input masked stack after geometry corrections and range
%       -looks = looks indexation 
%       -range_index = range vector indexation
%       -nf_p = non-fitting parameters strucrture
%       -cnf_p_ACDC = configuration parameters for fitting
%       -SWH = estimated SWH from the preliminary fitting
%       -epoch = estimated epoch from the preliminary fitting
%      OPTIONAL
%       -func_f0 = tabulated values of the func_f0
%       -func_f1 = tabulated values of the func_f1
%       
% OUTPUT:
%       -stack_ac = stack with the amplitude compensation performed
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - Depending on the stack being analyzed (first surfaces) is not always
% possible to have all the left and right beams: truncate to the available
% possible ones avoiding those set
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0:
% V1.1: 

function [stack_ac,k,l,g]=amplitude_compensation(stack,looks,range_index,nf_p,cnf_p_ACDC,SWH,epoch)


global N_samples zp_fact_range_cnf

%% ------------------- VARIABLES INITIALIZATION ---------------------------
%------------------ matrix notation of look and range indexes -------------
l=looks*ones(1,N_samples*zp_fact_range_cnf); %looks matrix
x=ones(nf_p.Neff,1)*range_index; %range matrix

%------------------ MSS initialization ------------------------------------
rou=cnf_p_ACDC.ini_rou;

%------------dilation term-------------------------------------------------
sigmaz  =   SWH/4.0; 
g       =   sqrt(2*nf_p.alphag_a*nf_p.alphag_r./(nf_p.alphag_a + nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * l.^2 + 2 * nf_p.alphag_a * nf_p.alphag_r * (sigmaz/nf_p.Lz)^2 ));        
    
%--------------Varibales  -------------------------------------------------
k=(x - epoch); xl=nf_p.Lx.*l; alpha_sigma=1/(nf_p.h^2*rou); k_prima=k; yk=nf_p.Ly.*abs(sqrt(k_prima));
            
%% --------- COMPUTATION OF THE ANTENNA/SURFACE ---------------------------
%Constant Term
Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2-alpha_sigma * xl.^2-nf_p.alphay * nf_p.yp^2-(nf_p.alphay + alpha_sigma).*(yk).^2).*cosh(2*nf_p.alphay*nf_p.yp*yk);
%Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2);

%% -------------- COMPENSATION OF BKL AND FUNF0 ---------------------------
stack_ac=stack./(sqrt(g).*Bkl);


end


