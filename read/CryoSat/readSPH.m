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
% READSPH: function that reads out the Specific Product Header of the specified file name
%
% Calling
%   [ sph ] = readSPH( filename )
%
% ----------------------------------------------------------
% 
% 
% Author:   Daniel Martinez / Pildo Labs
%           Josep Montolio / Pildo Labs
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Josep Montolio / Pildo Labs (20/10/09)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ sph ] = readSPH( filename )

fid = fopen(filename,'r','b');
mph = readMPH( fid );

size = mph.sph_size;
num_dsd = mph.num_dsd;

raw = fread(fid,size,'uchar');

p = 1;
%Prod description & identification
p = p +15;
p = p +1 ;
%added by EM 07.03.2016
product_id_str=char(raw(p:p+27)');
p = p +28;
p = p +1 ;
p = p +1 ;

%Prod Time info
p = p +22;
p = p +1 ;
%added by EM 07.03.2016
product_start_time=char(raw(p:p+26)');
p = p +27;
p = p +1 ;
p = p +1 ;
p = p +21;
p = p +1 ;
%added by EM 07.03.2016
product_stop_time=char(raw(p:p+26)');
p = p +27;
p = p +1 ;
p = p +1 ;
%added by EM 07.03.2016
product_time_info=struct('produc_start_time',product_start_time,...
                        'produc_stop_time',product_stop_time);

product_info=struct('product_id',char(mph.product'),...
                    'product_id_sph',product_id_str,...
                    'product_time_info',product_time_info);                    
%Prod Orbit Info
p = p +16;
% orbital phase code information
ABS_Orbit_Start=(char(raw(p:p+5)')); %added by EM 04.03.2016
p = p +6 ;
p = p +1 ;
p = p +24;
Rel_Time_ASC_Node_Start=(char(raw(p:p+10)')); %added by EM 04.03.2016
p = p +11;
p = p +3 ;
p = p +1 ;
p = p +15;
ABS_Orbit_Stop=(char(raw(p:p+5)')); %added by EM 04.03.2016
p = p +6 ;
p = p +1 ;
p = p +23;
Rel_Time_ASC_Node_Stop=(char(raw(p:p+10)')); %added by EM 04.03.2016
p = p +11;
p = p +3 ;
p = p +1 ;
p = p +23;
p = p +1 ;
Equator_Cross_Time=char(raw(p:p+26)'); %added by EM 04.03.2016
p = p +27;
p = p +1 ;
p = p +1 ;
p = p +19;
Equator_Cross_Long=(char(raw(p:p+10)')); %added by EM 04.03.2016
p = p +11;
p = p +10;
p = p +1 ;
p = p +15;
Ascending_Flag=char(raw(p)'); %added by EM 04.03.2016
p = p +1 ;
p = p +1 ;
% modified by EM 04.03.2016
%include the orbital info extracted from SPH
orbit_info=struct('phase_code',mph.orbit.phase_code,...
                  'cycle_num',mph.orbit.cycle_num,...
                  'rel_orbit',mph.orbit.rel_orbit,...
                  'ABS_Orbit_Start',ABS_Orbit_Start,...
                  'Rel_Time_ASC_Node_Start',Rel_Time_ASC_Node_Start,...
                  'ABS_Orbit_Stop',ABS_Orbit_Stop,...
                  'Rel_Time_ASC_Node_Stop',Rel_Time_ASC_Node_Stop,...
                  'Equator_Cross_Time',Equator_Cross_Time,...
                  'Equator_Cross_Long',Equator_Cross_Long,...
                  'Ascending_Flag',Ascending_Flag);


%Prod Location Info
p = p +10;
p = p +11;
p = p +10;
p = p +1 ;
p = p +11;
p = p +11;
p = p +10;
p = p +1 ;
p = p +9 ;
p = p +11;
p = p +10;
p = p +1 ;
p = p +10;
p = p +11;
p = p +10;
p = p +1 ;
p = p +50;
p = p +1 ;

%L0 Quality Info
p = p +13;
p = p +1 ;
p = p +1 ;
p = p +22;
p = p +6 ;
p = p +7 ;
p = p +1 ;
p = p +15;
p = p +6 ;
p = p +7 ;
p = p +1 ;
p = p +13;
p = p +1 ;
p = p +1 ;
p = p +12;
p = p +8 ;
p = p +1 ;
% cryosat2
%p = p +50;
%p = p +1 ;

p = p + 37;
p = p + 1;

%SIRAL Instr Config
p = p + 9;
p = p + 1;
ins_id = char(raw(p));%intrument identifier  A=Siral Nominal  B=Siral Redundant
p = p + 1;
p = p + 1;
p = p + 1;
p = p +12;
p = p +1 ;
p = p +10;
p = p +1 ;
p = p +1 ;
p = p +18;
p = p +1 ;
%added by EM 07.03.2016
ins_conf=char(raw(p:p+6)');
ins_info=struct('ins_id',ins_id,...
                'ins_conf',ins_conf);

% gnss & Doris sensor
doris_info=[];
gnss_info=[];

p = p +7 ;
p = p +1 ;
p = p +1 ;

% HACK: the L1 specification does not say it, but HERE there are more
% paremeters. We skip them
p = p + 179;

%L1 Proc Info
p = p +16;
p = p +1 ;
p = p +1 ;
p = p +14;
p = p +1 ;
p = p +1 ;
p = p +23;
p = p +6 ;
p = p +7 ;
p = p +1 ;
p = p +16;
p = p +6 ;
p = p +7 ;
p = p +1 ;
p = p +50;
p = p +1 ;


% DSD Section
%pre-allocation
dsds = repmat(struct(   'ds_name', '','ds_type','','filename', '','ds_offset',0,'ds_size',0,'num_dsr',0,'dsr_size',0),[num_dsd,1]);


for i = 1:num_dsd
	p = p +8 ;
	p = p +1 ;
    ds_name = char(raw(p:p + 27)');
	p = p +28;
	p = p +1 ;
	p = p +1 ;
	p = p +8 ;
    ds_type = char(raw(p));
	p = p +1 ;
	p = p +1 ;
	p = p +9 ;
	p = p +1 ;
    filename = char(raw(p:p + 61)');
	p = p +62;
	p = p +1 ;
	p = p +1 ;
	p = p +10;
    ds_offset = str2double(char(raw(p:p + 20)'));
	p = p +21;
	p = p +7 ;
	p = p +1 ;
	p = p +8 ;
    ds_size = str2double(char(raw(p:p + 20)'));
	p = p +21;
	p = p +7 ;
	p = p +1 ;
	p = p +8 ;
    num_dsr = str2double(char(raw(p:p + 10)'));
	p = p +11;
	p = p +1 ;
	p = p +9 ;
    dsr_size = str2double(char(raw(p:p + 10)'));
	p = p +11;
	p = p +7 ;
	p = p +1 ;
	p = p +32;
	p = p +1;
    dsds(i) = struct(   'ds_name', ds_name, ...
                        'ds_type',ds_type,...
                        'filename', filename,...
                        'ds_offset',ds_offset,...
                        'ds_size',ds_size,...
                        'num_dsr',num_dsr,...
                        'dsr_size',dsr_size);
end




%added EM 04.03.2016 & Modified EM 07.03.2016
sph = struct(   'raw', raw,...
                'dsds', dsds,...
                'ins_info',ins_info,...
                'doris_info',doris_info,...
                'gnss_info',gnss_info,...
                'orbit_info',orbit_info,...
                'product_info',product_info,...
                'acq_station',mph.acq_station); 

% disp( char(sph.raw') );