% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% CR2 SARIn 
% ---------------------------------------------------------
% Objective: Define Constants parameters
% 
% INPUTs : - 
% OUTPUTs: -
%
% ----------------------------------------------------------
% Author:    Roger Escolà  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Mònica Roca   / isardSAT
% Last rev.: Mònica Roca   / isardSAT (11/09/2013)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global c_cst pi_cst sec_in_day_cst mean_sat_alt_cst
global earth_radius_cst flat_coeff_cst semi_major_axis_cst semi_minor_axis_cst

semi_minor_axis_cst = 6356752.3142;
mean_sat_alt_cst = 717242;
sec_in_day_cst = 86400;
c_cst = 299792458;
pi_cst = pi;
flat_coeff_cst = 0.00335281066;
earth_radius_cst = 6378137;
semi_major_axis_cst = 6378137;

disp('Inside CST')
