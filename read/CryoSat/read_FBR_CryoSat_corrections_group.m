function [burst] = read_FBR_CryoSat_corrections_group(fid, burst)


        %-----------------------------%
        %--  read Corrections Group --%
        %-----------------------------%
        %35 Dry Tropospheric Correction mm 4 sl
        burst.dry_tropo_correction_bursts = fread(fid,1,'int32')  * 1e-3;
        %36 Wet Tropospheric Correction mm 4 sl
        burst.wet_tropo_correction_bursts = fread(fid,1,'int32') * 1e-3;
        %37 Inverse Barometric Correction mm 4 sl
        burst.inverse_baro_correction_bursts = fread(fid,1,'int32') * 1e-3;
        %38 Dynamic Atmospheric Correction mm 4 sl
        burst.Dynamic_atmospheric_correction_bursts = fread(fid,1,'int32') * 1e-3;        
        %39 GIM Ionospheric Correction mm 4 sl
        burst.GIM_iono_correction_bursts = fread(fid,1,'int32') * 1e-3;
        %40 Model Ionospheric Correction mm 4 sl
        burst.model_iono_correction_bursts = fread(fid,1,'int32') * 1e-3;
        %41 Ocean Equilibrium Tide mm 4 sl
        burst.ocean_equilibrium_tide_bursts = fread(fid,1,'int32') * 1e-3;
        %42 Long Period Tide Height 4 sl
        burst.long_period_tide_height_bursts = fread(fid,1,'int32') * 1e-3;
        %43 Ocean Loading Tide mm 4 sl
        burst.ocean_loading_tide_bursts = fread(fid,1,'int32') * 1e-3;
        %44 Solid Earth Tide mm 4 sl
        burst.solid_earth_tide_bursts = fread(fid,1,'int32') * 1e-3;
        %45 Geocentric Polar Tide mm 4 sl
        burst.geocentric_polar_tide_bursts = fread(fid,1,'int32') * 1e-3;
        %46 Surface type flag 4 ul
        burst.surface_type_flag_bursts = fread(fid,1,'uint32');
        fseek(fid,4+4+4+4,'cof');
%         %47 spares 4*1 uc
%         fread(fid,4,'uint8');
%         %48 Correction status flags 4 ul (see table 2.3.3-5)
%         %correction_status_flags = fread(fid,1,'uint32');
% 
%         correction_status_flags(1) = fread(fid,1,'ubit1');
%         correction_status_flags(2) = fread(fid,1,'ubit1');
%         correction_status_flags(3) = fread(fid,1,'ubit1');
%         correction_status_flags(4) = fread(fid,1,'ubit1');
%         correction_status_flags(5) = fread(fid,1,'ubit1');
%         correction_status_flags(6) = fread(fid,1,'ubit1');
%         correction_status_flags(7) = fread(fid,1,'ubit1');
%         correction_status_flags(8) = fread(fid,1,'ubit1');
%         correction_status_flags(9) = fread(fid,1,'ubit1');
%         correction_status_flags(10) = fread(fid,1,'ubit1');
%         correction_status_flags(11) = fread(fid,1,'ubit1');
%         correction_status_flags(12) = fread(fid,1,'ubit1');
%         correction_status_flags(13) = fread(fid,1,'ubit20');
%                 
%         
%         %49 Correction error flags 4 ul (see table 2.3.3-6)
%         %correction_error_flags = fread(fid,1,'uint32');
%         
%         correction_error_flags(1) = fread(fid,1,'ubit1');
%         correction_error_flags(2) = fread(fid,1,'ubit1');
%         correction_error_flags(3) = fread(fid,1,'ubit1');
%         correction_error_flags(4) = fread(fid,1,'ubit1');
%         correction_error_flags(5) = fread(fid,1,'ubit1');
%         correction_error_flags(6) = fread(fid,1,'ubit1');
%         correction_error_flags(7) = fread(fid,1,'ubit1');
%         correction_error_flags(8) = fread(fid,1,'ubit1');
%         correction_error_flags(9) = fread(fid,1,'ubit1');
%         correction_error_flags(10) = fread(fid,1,'ubit1');
%         correction_error_flags(11) = fread(fid,1,'ubit1');
%         correction_error_flags(12) = fread(fid,1,'ubit1');
%         correction_error_flags(13) = fread(fid,1,'ubit20');
%         %50 Spare 4*1 uc
%         fread(fid,4,'uint8');
end