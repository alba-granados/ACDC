% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT Ltd. 
% --------------------------------------------------------
%
% Ground Processor Prototype for altimetry missions:
% CryoSat-2 SARIn
% CryoSat-2 SAR
% Sentinel 6 RAW
% Sentinel 6 RMC
% 
% ---------------------------------------------------------
% Inputs: 
% L0 (ISP)+ orbit file + attitude file
% L1A (FBR in case of CryoSat-2
% L1B-S
% L1B
% 
% Output
%   L1A record
%
% ----------------------------------------------------------
% 
% Author:   Albert Garcia-Mondejar / isardSAT
%           Roger Escola Jane / isardSAT
%
% Reviewer: Monica Roca / isardSAT
%
% Last revision: Monica Roca / isardSAT (13/05/11)
%
% This software is subject to the conditions 
% set forth in the ESA contract XXXXXXXXXX
% v1.1 The file identifier is passed instead the name to avoid opening the file 
% each time a burst is read
% v2.0 2016/05/13 Sentinel 3 added. Deleted reading of SPH in CR2. fid is copied in the files struct
% v2.1 2016/05/18 Removed ds_offset and ds_size as globals
% v2.2 2016/11/16 fclose moved from L1B_processing chain to here, N_burst added as input to check when we read the last burst
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L1A,files] = read_L1A_record(files,L1A,original_burst, i_burst,N_bursts)

global mode  mission


switch mission
    case 'S3_'
        if i_burst==1
            files.fid = netcdf.open(files.filename_L1A,'NC_NOWRITE'); %open file
        end
        [netCDF_L1A] = readanyNETCDF_record(files.fid,original_burst);
        [L1A]        = adapt_netCDF2internal(netCDF_L1A,files.filename_L1A,original_burst);
        if i_burst==N_bursts
            netcdf.close(files.fid); %close the netcdf file
        end
    case 'CR2'
        if i_burst==1

            % Open the binary file only once and then use the pointer to avoid
            % several open/close iterations per burst
            files.fid = fopen(files.filename_L1A,'r','b');
        end
        
        [L1A] = readFBR_CryoSat_SAR_record(files,L1A, original_burst);
        
        if i_burst==N_bursts
            fclose(files.fid); %close the input binary file
        end
 end
    %fclose all;
end

