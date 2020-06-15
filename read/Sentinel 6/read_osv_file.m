 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   READING Jason_CS OSV   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   isardSAT S.L.         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coordinate_system,datation_tai,...
    position_x,position_y,position_z,position_modulus,...
    velocity_x,velocity_y,velocity_z,velocity_modulus] = read_osv_file(Pathin)

     
   
  

% tic

% L1bfolder=dir(JC_IC4D_FSV_001.dat);

% for i=3:size(L1bfolder,1)
% 
% filename=L1bfolder(i,1).name

% -- jump headers --
elements_size=20;
%sphdsd=6975;

elements_offset=201;
mph=1400;


% read the number of dsr in the dsd 
fid = fopen(Pathin,'r','b');
    fseek(fid, elements_offset, 'bof');
%     ftell(fid)

    numdsr_strchar=fread(fid,elements_size,'*char');
    numdsr_str=[];
    %numdsr=cellstr(numdsr_strchar);
   for i_index = 1:length(numdsr_strchar)
     numdsr_str = [numdsr_str,numdsr_strchar(i_index)];
   end
     
         
   numdsr=str2num(numdsr_str);
     
   fseek(fid,mph, 'bof');    
    
  datation_tai=zeros(numdsr,1);
     position_x=zeros(numdsr,1);
     position_y=zeros(numdsr,1);
     position_z=zeros(numdsr,1);
     position_modulus=zeros(numdsr,1);
     velocity_x=zeros(numdsr,1);
     velocity_y=zeros(numdsr,1);
     velocity_z=zeros(numdsr,1);
     velocity_modulus=zeros(numdsr,1);
     velocity_modulus=zeros(numdsr,1);
          
for i=1:numdsr
    
    %----------------------------%
    %--      read science      --%
    %----------------------------%
  
	
	%1 Coordinate system 
         
 coordinate_system(i) = fread(fid,1,'int32');
%           
%    %2 Derivatives
% derivatives(i) = fread(fid,1,'int32');
%    %3 Distance Unit
% distance_unit(i)= fread(fid,3,'int32');
%     4 Reliability
% reliabbility(i) = fread(fid,1,'int32');

        fseek(fid,12, 'cof');
%     %5 datation TAI
        datation_tai(i) = fread(fid,1,'double' );
%     %6 Transport datation TAI 1
% transport_datation_tai_1(i) = fread(fid,1,'double' );
      %7 Transport datation TAI 2
% transport_datation_tai_1(i) = fread(fid,1,'double' );
      %8 Transport datation TAI 3
% transport_datation_tai_1(i) = fread(fid,1,'double' );
      %9 Transport datation TAI 4
% transport_datation_tai_1(i) = fread(fid,1,'double' );

        fseek(fid,32, 'cof');
 
 %     %10 Position X
position_x(i) = fread(fid,1,'double' );
 %     %11 Position Y
position_y(i) = fread(fid,1,'double' );
 %     %12 Position Z
position_z(i) = fread(fid,1,'double' );
 %     %13 Position modulus
position_modulus(i) = fread(fid,1,'double' );
 %     %14 Velocity X
velocity_x(i) = fread(fid,1,'double' );
 %     %15 Velocity Y
velocity_y(i) = fread(fid,1,'double' );
 %     %16 Velocity Z
velocity_z(i) = fread(fid,1,'double' );
 %     %17 Velocity modulus
velocity_modulus(i) = fread(fid,1,'double' );
 %     %18 Spares
% spare_1(i) = fread(fid,1,'d' );
 %     %19 Spares
% spare_2(i) = fread(fid,1,'d' );
 %     %20 Spares
% spare_3(i) = fread(fid,1,'d' );
 %     %21 Spares
% spare_4(i) = fread(fid,1,'d' );

    

fseek(fid, 32, 'cof');

end

    fclose(fid);
%     toc
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%      ----------------- SCIENCE DATA ---------------------	Byte length
% Identifier                    Description                     Unit                Type	Size
%1                              Coordinate system=4         	N/A                 sl      4
%2                              Derivatives=1                   N/A                 sl      4
%3                          	Distance unit=0                 N/A                 sl      4
%4                              Reliability=1                   N/A                 sl      4
%5                              Datation (TAI)                  s                   d       8
%6                              Transport datation 1[1]         TBD                 d       8
%7                              Transport datation 2                                d       8
%8                              Transport datation 3                                d       8
%9                              Transport datation 4                                d       8
%10                             Position x                      m                   d       8
%11                             Position y                      m                   d       8
%12                             Position z                      m                   d       8
%13                             Position modulus                m                   d       8
%14                             Velocity x                      m/s                 d       8
%15                             Velocity y                      m/s                 d       8
%16                             Velocity z                      m/s                 d       8
%17                             Velocity modulus                m/s                 d       8
%18                             Spare                           N/A                 d       8
%19                             Spare                           N/A                 d       8
%20                             Spare                           N/A                 d       8
%21                             Spare                           N/A                 d       8
%				Total 152
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%