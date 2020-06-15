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
% MODULE(vec) 
%   Calculates the module of _vec_
%
% MODULE(mat)
%   Calculates the module of the rows
%   of _mat_	
%
% A.J Quiles, 3-8-97
% $ Revision 1.0, 1-10-98
% $ 	also takes matrix as arguments
% $ Revision 1.1, 05-10-98
% $	admits row vectors
% $ Revision 2.0, 24.02.05 - uses rows in matrixes
% 
% Author:   Josep Montolio / Pildo Labs
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Josep Montolio / Pildo Labs (19/10/09)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mod = module(a)


if nargin < 1
	help module
	error('Missing Arguments.')
end

% Checks wether it is a column vector
if size(a, 2) == 1, a = a'; end

for nv = 1:size(a, 1)
	vec = a(nv, :);
	mod(nv,1) = sqrt(vec*vec');
end