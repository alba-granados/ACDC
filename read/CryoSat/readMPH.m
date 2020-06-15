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
% READMPH: function that reads out the Main Product Header from current
% position on the file_handler
%
% Calling
%   [ mph ] = readMPH( file_handler )
%
% ----------------------------------------------------------
% 
% 
% Author:   Mercedes Reche / Pildo Labs
%           Josep Montolio / Pildo Labs
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Josep Montolio / Pildo Labs (20/10/09)
%
% $Id: readMPH.m 152 2005-12-20 09:29:33Z  $
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ mph ] = readMPH( file_handler )

raw = fread(file_handler, 1247,'uchar');

p = 1; 
%PRODUCT field
%PRODUCT identification Info (comment added by EM 19.02.2016)
p = p + 8  ;
p = p + 1  ;
product = raw(p:p + 61);
p = p + 62 ;

p = p + 1  ;
p = p + 1  ;
p = p + 11 ;
p = p + 1  ;
p = p + 1  ;
p = p + 8  ;
p = p + 1  ;
p = p + 23 ;
p = p + 1  ;
p = p + 1  ;
p = p + 40 ;
p = p + 1  ;
% Data Processing Information (comment added by EM 19.02.206)
p = p + 20 ;
p = p + 1  ;
acq_station=deblank(strvcat(raw(p:p+19)'));
p = p + 20 ;
p = p + 1  ;
p = p + 1  ;
p = p + 12 ;
p = p + 1  ;
p = p + 6  ;
p = p + 1  ;
p = p + 1  ;
p = p + 10 ;
p = p + 1  ;
p = p + 27 ;
p = p + 1  ;
p = p + 1  ;
p = p + 13 ;
p = p + 1  ;
p = p + 14 ;
p = p + 1  ;
p = p + 1  ;
p = p + 40 ;
p = p + 1  ;
%Information on Time of Data (comment added by EM 19.02.2016)
p = p + 14 ;
p = p + 1  ;
p = p + 27 ;
p = p + 1  ;
p = p + 1  ;
p = p + 13 ;
p = p + 1  ;
p = p + 27 ;
p = p + 1  ;
p = p + 1  ;
p = p + 40 ;
p = p + 1  ;
%Orbit Information (comment added by EM 19.02.2016)
p = p + 6  ;
% orbital phase code information
phase_code=(char(raw(p)')); %added by EM 19.02.2016
p = p + 1  ;
p = p + 1  ;
p = p + 6  ;
% orbital cycle information
cycle_num=char(raw(p:p+3)'); %added by EM 19.02.2016
p = p + 4  ;
p = p + 1  ;
p = p + 10 ;
% relative orbit number
rel_orbit=char(raw(p:p+5)'); %added by EM 19.02.2016
p = p + 6  ;
p = p + 1  ;
p = p + 10 ;
% absolute orbit number
abs_orbit=char(raw(p:p+5)'); %added by EM 19.02.2016
%Generate the struct info with the orbital information required by Upporto
% added by EM 19.02
orbit=struct('phase_code',phase_code,...
                  'cycle_num',cycle_num,...
                  'rel_orbit',rel_orbit,...
                  'abs_orbit',abs_orbit);
p = p + 6  ;
p = p + 1  ;
p = p + 18 ;
p = p + 1  ;
p = p + 27 ;
p = p + 1  ;
p = p + 1  ;
p = p + 10 ;
p = p + 8  ;
p = p + 3  ;
p = p + 1  ;
p = p + 11 ;
p = p + 12 ;
p = p + 3  ;
p = p + 1  ;
p = p + 11 ;
p = p + 12 ;
p = p + 3  ;
p = p + 1  ;
p = p + 11 ;
p = p + 12 ;
p = p + 3  ;
p = p + 1  ;
p = p + 11 ;
p = p + 12 ;
p = p + 5  ;
p = p + 1  ;
p = p + 11 ;
p = p + 12 ;
p = p + 5  ;
p = p + 1  ;
p = p + 11 ;
p = p + 12 ;
p = p + 5  ;
p = p + 1  ;
p = p + 14 ;
p = p + 1  ;
p = p + 2  ;
p = p + 1  ;
p = p + 1  ;
p = p + 40 ;
p = p + 1  ;
p = p + 13 ;
p = p + 1  ;
p = p + 27 ;
p = p + 1  ;
p = p + 1  ;
p = p + 16 ;
p = p + 11 ;
p = p + 1  ;
p = p + 11 ;
p = p + 11 ;
p = p + 4  ;
p = p + 1  ;
p = p + 32 ;
p = p + 1  ;
p = p + 9  ;
p = p + 1  ;
p = p + 27 ;
p = p + 1  ;
p = p + 1  ;
p = p + 10 ;
p = p + 4  ;
p = p + 1  ;
p = p + 9  ;
p = p + 1  ;
p = p + 1  ;
p = p + 40 ;
p = p + 1  ;
p = p + 12 ;
p = p + 1  ;
p = p + 1  ;
p = p + 9  ;
p = p + 21 ;
p = p + 7  ;
p = p + 1  ;
p = p + 9  ;

%SPH_SIZE
sph_size = str2double(char(raw(p:p + 10)'));

p = p + 11 ;
p = p + 7  ;
p = p + 1  ;
p = p + 8  ;

%NUM_DSD
num_dsd = str2double(char(raw(p:p + 10)'));

p = p + 11 ;
p = p + 1  ;
p = p + 9  ;
p = p + 11 ;
p = p + 7  ;
p = p + 1  ;
p = p + 14 ;
p = p + 11 ;
p = p + 1  ;
p = p + 4  ;
p = p + 6  ;
p = p + 1  ;
p = p + 29 ;
p = p + 1  ;


mph =   struct( 'raw',raw,...
                'product', product,...
                'acq_station',acq_station,...
                'sph_size',sph_size,...
                'num_dsd',num_dsd,...
                'orbit',orbit); % modified by EM 19.02.2016 to include the orbital information required UPorto
