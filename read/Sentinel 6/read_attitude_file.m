%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   READING Jason_CS Attitude   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   isardSAT S.L.         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [datation_tai, roll, pitch, yaw] = read_attitude_file(Pathin)

 
global roll_random_error_chd pitch_random_error_chd yaw_random_error_chd
global roll_bias_error_chd pitch_bias_error_chd yaw_bias_error_chd
global roll_harmonic_error_chd pitch_harmonic_error_chd yaw_harmonic_error_chd
global random_pointing_error_cnf bias_pointing_error_cnf harmonic_pointing_error_cnf


% tic

% -- jump headers --
elements_size = 20;
elements_offset = 201;
mph = 1400;


% read the number of dsr in the dsd
fid = fopen(Pathin,'r','b');
fseek(fid, elements_offset,'bof');
% ftell(fid)

numdsr_strchar = fread(fid,elements_size,'*char');
numdsr_str = [];

for i_index = 1:length(numdsr_strchar)
	numdsr_str = [numdsr_str,numdsr_strchar(i_index)];
end
numdsr = str2double(numdsr_str);

fseek(fid,mph,'bof');

datation_tai = zeros(numdsr,1);
roll = zeros(numdsr,1);
pitch = zeros(numdsr,1);
yaw = zeros(numdsr,1);

for i=1:numdsr
    
    %----------------------------%
    %--      read science      --%
    %----------------------------%
    
	
	%1 Datation TAI
    datation_tai(i) = fread(fid,1,'double');
    
    fseek(fid,96,'cof');
    
    %14 Roll
    roll(i) = fread(fid,1,'double');
    %15 Pitch
    pitch(i) = fread(fid,1,'double');
    %16 Yaw
    yaw(i) = fread(fid,1,'double');
    
end

%% ADDING RANDOM ERROR;

% roll    = roll  + random_pointing_error_cnf(1)*roll_random_error_chd*randn(numdsr,1)+ ...
%                   bias_pointing_error_cnf(1)*roll_bias_error_chd +...
%                   harmonic_pointing_error_cnf(1)*roll_harmonic_error_chd;
% pitch   = pitch + random_pointing_error_cnf(2)*pitch_random_error_chd*randn(numdsr,1)+ ...
%                   bias_pointing_error_cnf(2)*pitch_bias_error_chd+...
%                   harmonic_pointing_error_cnf(2)*pitch_harmonic_error_chd;
% yaw     = yaw   + random_pointing_error_cnf(3)*yaw_random_error_chd*randn(numdsr,1)+ ...
%                   bias_pointing_error_cnf(3)*yaw_bias_error_chd+...
%                   harmonic_pointing_error_cnf(3)*yaw_harmonic_error_chd;

fclose(fid);


end
    