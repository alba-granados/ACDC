% READ METEO DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% S3 MPC 
% This code implements the READING for the meteo rinex format
% INPUTs: 
%
%     Name                      Units	Type	Origin
%     Filename         
%     Rinex file
%       ftp://igs.org/pub/data/format/rinex302.pdf
%
%       # / TYPES OF OBSERV - Number of different observation types stored in the file
%       - Observation types - The following meteorological observation types are defined in RINEX Version 3:
%       - PR : Pressure (mbar)
%       - TD : Dry temperature (deg Celsius)
%       - HR : Relative humidity (percent)
%       - ZW : Wet zenith path delay (mm), (for WVR data)
%       - ZD : Dry component of zen.path delay (mm)
%       - ZT : Total zenith path delay (mm)
%       - WD : Wind azimuth (deg) from where the wind blows
%       - WS : Wind speed (m/s)
%       - RI : "Rain increment" (1/10 mm): Rain accumulation since last measurement
%       - HI : Hail indicator non-zero: Hail detected since last measurement
%       - The sequence of the types in this record must correspond to the sequence of the measurements in the data records.
%       - If more than 9 observation types are being used, use continuation lines with format (6X,9(4X,A2)) 
% OUTPUTs:
%     Name                      Units	Type	Destination
%     matlab structure with the meteo data

function [meteo]=read_meteo_rinex(filename)

Row_offset  = 0;
meteo       = struct();
fid_meteo   = fopen(filename,'r');

while ~feof(fid_meteo)
    tline       = fgetl(fid_meteo);
    %create structure fields from columns of RINEX file
    Row_offset=Row_offset+1;
    if(strfind(tline,'TYPES OF OBSERV'))
        spaces  = strfind(tline,'    ');
        N_params = str2double(tline(spaces(3)-1));
        for i_param = 1:N_params
            param(i_param,:)= tline(spaces(3+i_param)-2:spaces(3+i_param)-1);
            
        end
        
    end    
    if(strfind(tline,'END OF HEADER'))
        break;        
    end
    
end
fclose(fid_meteo);
aux = dlmread(filename,'',Row_offset,1);

meteo.day       = aux(:,1);
meteo.month     = aux(:,2);
meteo.year      = aux(:,3);
meteo.hour      = aux(:,4);
meteo.second    = aux(:,5);
for i_param = 1:N_params
    meteo.(param(i_param,:))=aux(:,5+i_param);
end


end

