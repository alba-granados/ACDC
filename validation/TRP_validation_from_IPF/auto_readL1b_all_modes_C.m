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
% READL1b: function that reads the Level1b data set from input filename. 
%
% Calling
%   out_L1B = auto_readL1b_all_modes (filename, headers, KML)
%
% Inputs
%   filename: input LRM, SAR or SARIN L1b file
%   headers:   if headers = false --> reads num_records from file assuming no headers
%              if headers = true --> reads all records specified in file headers (sph)
%   num_records: number of records to be read
%
% Output
%   l1b_ds           : data contained in the file 
%
% Comments:
%   - num_records input is not used due the function has been automated in
%     order to take all the available records.
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
%           Daniel Martínez / Pildo Labs
%           Josep Montolio / Pildo Labs
%           Pablo García / isardSAT
%           Albert García / isardSAT
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Mònica Roca / isardSAT (26/05/11)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [l1b_ds,num_dsr] = auto_readL1b_all_modes_C( filename, headers, KML)



if ( headers )
%     disp(['reading L1b from ' filename ' WITH headers'] );
    sph = readSPH(filename);

    % L1b data set must be the first one
    l1b_dsd = sph.dsds(1);
    % Siral Identificator
   
    %Looking for Mode
    if(strfind(l1b_dsd.ds_name,'CAL1_LRM'))
        disp('SIR_L1B_CAL1_LRM data set found');
    elseif(strfind(l1b_dsd.ds_name,'CAL1_SARIN'))
        disp('SIR_L1B_CAL1_SARIN data set found');
    elseif(strfind(l1b_dsd.ds_name,'CAL1_SAR'))
        disp('SIR_L1B_CAL1_SAR data set found');
    elseif(strfind(l1b_dsd.ds_name,'CAL2_SARIN'))
        disp('SIR_L1B_CAL2_SARIN data set found');
    elseif(strfind(l1b_dsd.ds_name,'CAL2_SAR'))
        disp('SIR_L1B_CAL2_SAR data set found');
    elseif(strfind(l1b_dsd.ds_name,'TRK_CAL3'))
        disp('SIR_L1B_TRK_CAL3 data set found');
    elseif(strfind(l1b_dsd.ds_name,'TRK_SARIN'))
        disp('SIR_L1B_TRK_SARIN data set found');
    elseif(strfind(l1b_dsd.ds_name,'TRK_SAR'))
        disp('SIR_L1B_TRK_SAR data set found');
    elseif(strfind(l1b_dsd.ds_name,'CAL3'))
        disp('SIR_L1B_CAL3 data set found');
    elseif(strfind(l1b_dsd.ds_name,'CAL4'))
        disp('SIR_L1B_CAL4 data set found');
    elseif(strfind(l1b_dsd.ds_name,'ACQ'))
        disp('SIR_L1B_ACQ data set found');
    elseif(strfind(l1b_dsd.ds_name,'LRM'))
        siralMode= 'LRM';
        disp('SIR_L1B_LRM data set found');
    elseif(strfind(l1b_dsd.ds_name,'SARIN'))
        siralMode= 'SARin';
%         disp('SIR_L1B_SARIN data set found');
    elseif(strfind(l1b_dsd.ds_name,'SAR'))
        siralMode= 'SAR';
        disp('SIR_L1B_SAR data set found');
    
    else
        disp('not valid data found');
        return;
    end

    ds_offset = l1b_dsd.ds_offset;
    ds_size = l1b_dsd.ds_size;
    num_dsr = l1b_dsd.num_dsr;
    dsr_size = l1b_dsd.dsr_size;
else
    disp(['reading L1b from ' filename ' WITHOUT headers'] );
   
    ds_offset = 0;%In case there is no header
end

% read the dsrs
fid = fopen(filename,'r','b');
fseek(fid, ds_offset, 0);
%num_dsr=100;
%iterate available dsr's

for iRecord=1:num_dsr
    
    % Progress bar
%     customText = 'Loop...';
%     percentageDone =  iRecord / num_dsr;
%     stopBar= progressbar(percentageDone, 0, customText);
%     if (stopBar) 
%         break; 
%     end
    
    %% read Time and Orbit Group %
   
    
    for j=1:20 %20 items per record
              
		%1 Data Record Time (MDSR TimeStamp) TAI 12 sl+2*ul
        day = fread(fid,1,'int32');
        
%         if(isempty(day)==1)%If we don't find anything inside means that the file has end
%             end_of_file=1;
%             record_num = record_num-1;
%         else
                       
        second = fread(fid,1,'uint32');
        microsecond = fread(fid, 1, 'uint32');
        data_record_time = struct(  'day', day,...
                                    'second', second, ...
                                    'microsecond',microsecond);

        %2 USO Correction Factor (ratio) 10-9 4 sl
        USO_correction = fread(fid,1,'uint32');
        %3 Mode ID 2 us (see table 2.3.3-2) 
        modeID_mode = fread(fid,1,'ubit6');
        modeID_SARIn_deg_case = fread(fid,1,'ubit1');
        modeID_reserved = fread(fid,1,'ubit1');
        modeID_CAL4 = fread(fid,1,'ubit1');
        modeID_Platform_Cont = fread(fid,1,'ubit2');
        modeID_reserved2 = fread(fid,1,'ubit5'); 

        modeID =struct(    'modeID_mode', modeID_mode,...
                                'modeID_SARIn_deg_case', modeID_SARIn_deg_case,...
                                'modeID_reserved', modeID_reserved,...
                                'modeID_CAL4', modeID_CAL4,...
                                'modeID_Platform_Cont', modeID_Platform_Cont,...
                                'modeID_reserved2',modeID_reserved2);


        %4 Source Sequence Counter 2 us (see note 6)
        ssc = fread(fid,1,'uint16');
        %5 Instrument Configuration 4 ul (see table 2.3.3-3)

        ins_rx_in_use = fread(fid,1,'ubit2');
        ins_siral_id = fread(fid,1,'ubit1');
        ins_reserved = fread(fid,1,'ubit1');
        ins_bandwidth = fread(fid,1,'ubit2');
        ins_reserved2 = fread(fid,1,'ubit2');
        ins_tracking_mode = fread(fid,1,'ubit2');
        ins_ext_calibration = fread(fid,1,'ubit1');
        ins_reserved3 = fread(fid,1,'ubit1');
        ins_loop_stat = fread(fid,1,'ubit1');
        ins_loss_echo = fread(fid,1,'ubit1');
        ins_real_time_error = fread(fid,1,'ubit1');
        ins_echo_sat_error = fread(fid,1,'ubit1');
        ins_rx_band_attenuation = fread(fid,1,'ubit1');
        ins_cycle_report = fread(fid,1,'ubit1');
        ins_star_tracker_1 = fread(fid,1,'ubit1');
        ins_star_tracker_2 = fread(fid,1,'ubit1');
        ins_star_tracker_3 = fread(fid,1,'ubit1');
        ins_reserved4 = fread(fid,1,'ubit11');

        instrument_configuration =struct(    'ins_rx_in_use', ins_rx_in_use,...
                            'ins_siral_id', ins_siral_id,...
                            'ins_reserved', ins_reserved,...
                            'ins_bandwidth', ins_bandwidth,...
                            'ins_reserved2', ins_reserved2,...
                            'ins_tracking_mode',ins_tracking_mode,...                                
                            'ins_ext_calibration', ins_ext_calibration,...
                            'ins_reserved3', ins_reserved3,...
                            'ins_loop_stat', ins_loop_stat,...
                            'ins_loss_echo', ins_loss_echo,...
                            'ins_real_time_error',ins_real_time_error,...
                            'ins_echo_sat_error', ins_echo_sat_error,...
                            'ins_rx_band_attenuation', ins_rx_band_attenuation,...
                            'ins_cycle_report', ins_cycle_report,...
                            'ins_star_tracker_1', ins_star_tracker_1,...
                            'ins_star_tracker_2', ins_star_tracker_2,...
                            'ins_star_tracker_3', ins_star_tracker_3,...
                            'ins_reserved4',ins_reserved4);

        %6 Burst counter (always starts from 1 and incremented at group rate) 4 ul
        burst_counter = fread(fid,1,'uint32');        
        %7 Latitude of measurement 10-1 mcrodeg 4 sl (see note 1)
        latitude = fread(fid,1,'int32');  
        %8 Longitude of measurement 10-1 microdeg 4 sl (see note 1)
        longitude = fread(fid,1,'int32');  
        %9 Altitude of COG above reference ellipsoid (interpolated value)mm 4 sl
        altitude_COG = fread(fid,1,'int32');
        %10 Instantaneous altitude rate derived from orbit mm/s 4 sl
        instantaneous_altitude_rate = fread(fid,1,'int32'); 
        %11 Satellite velocity vector[3](in ITRF) mm/s 3*4 sl
        satellite_velocity(1) = fread(fid,1,'int32'); 
        satellite_velocity(2) = fread(fid,1,'int32');
        satellite_velocity(3) = fread(fid,1,'int32');
        %12 Real beam direction vector[3](in CRF) micros 3*4 sl
        real_beam_direction(1) = fread(fid,1,'int32');
        real_beam_direction(2) = fread(fid,1,'int32');
        real_beam_direction(3) = fread(fid,1,'int32');
        %13 Interferometer baseline vector[3](in CRF)micros 3*4 sl
        inferometer_baseline(1) = fread(fid,1,'int32');
        inferometer_baseline(2) = fread(fid,1,'int32');
        inferometer_baseline(3) = fread(fid,1,'int32');
        
        %14
        str_id      = fread(fid,1,'uint16');
        %15 Antenna Bench Roll Angle 10^-7 deg
        bench_roll  = fread(fid,1,'int32');
        %16 Antenna Bench Pitch Angle 10^-7 deg
        bench_pitch  = fread(fid,1,'int32');
        %17 Antenna Bench Yaw Angle 10^-7 deg
        bench_yaw   = fread(fid,1,'int32');
        
        
        
        %18 FBR Measurement ConfidenceData (flag word)4 u

        confi_block_degraded = fread(fid,1,'ubit1');
        confi_blank_block = fread(fid,1,'ubit1');
        confi_datation_degraded = fread(fid,1,'ubit1');
        confi_orbit_prop_error = fread(fid,1,'ubit1');
        confi_orbit_file_change = fread(fid,1,'ubit1');
        confi_orbit_discon = fread(fid,1,'ubit1');
        confi_echo_sat = fread(fid,1,'ubit1');
        confi_other_echo_error = fread(fid,1,'ubit1');
        confi_rx1_error_for_SARIN = fread(fid,1,'ubit1');
        confi_rx2_error_for_SARIN = fread(fid,1,'ubit1');
        confi_window_delay_inconsistency = fread(fid,1,'ubit1');
        confi_AGC_inconsistency = fread(fid,1,'ubit1');
        confi_cal1_correction_miss = fread(fid,1,'ubit1');
        confi_cal1_correction_from_IPF_DB = fread(fid,1,'ubit1');
        confi_DORIS_USO_correction = fread(fid,1,'ubit1');
        confi_complex_cal1_correction_from_IPF_DB = fread(fid,1,'ubit1');
        confi_TRK_echo_error = fread(fid,1,'ubit1');
        confi_echo_rx1_error = fread(fid,1,'ubit1');
        confi_echo_rx2_error = fread(fid,1,'ubit1');
        confi_NMP_inconsistency = fread(fid,1,'ubit1');
        confi_azimuth_cal_missing = fread(fid,1,'ubit1');
        confi_azimuth_cal_from_IPF_DB = fread(fid,1,'ubit1');
        confi_range_window_cal_function_missing = fread(fid,1,'ubit1');
        confi_range_window_cal_function_from_IPF_DB = fread(fid,1,'ubit1');
        confi_reserved = fread(fid,1,'ubit1');
        confi_cal2_correction_missing = fread(fid,1,'ubit1');
        confi_cal2_correction_from_IPF_DB = fread(fid,1,'ubit1');
        confi_power_scaling_error_LRM_only = fread(fid,1,'ubit1');
        confi_attitude_correction_missing = fread(fid,1,'ubit1');
        confi_reserved2 = fread(fid,1,'ubit2');
        confi_phase_perturbation_correction = fread(fid,1,'ubit1');
        spare  = fread(fid,1,'int32');
        L1B_measurement_confidence_flag =struct(    'confi_block_degraded', confi_block_degraded,...
                            'confi_blank_block', confi_blank_block,...
                            'confi_datation_degraded', confi_datation_degraded,...
                            'confi_orbit_prop_error', confi_orbit_prop_error,...
                            'confi_orbit_file_change', confi_orbit_file_change,...
                            'confi_orbit_discon',confi_orbit_discon,...                                
                            'confi_echo_sat', confi_echo_sat,...
                            'confi_other_echo_error', confi_other_echo_error,...
                            'confi_rx1_error_for_SARIN', confi_rx1_error_for_SARIN,...
                            'confi_rx2_error_for_SARIN', confi_rx2_error_for_SARIN,...
                            'confi_window_delay_inconsistency',confi_window_delay_inconsistency,...
                            'confi_AGC_inconsistency', confi_AGC_inconsistency,...
                            'confi_cal1_correction_miss', confi_cal1_correction_miss,...
                            'confi_cal1_correction_from_IPF_DB', confi_cal1_correction_from_IPF_DB,...
                            'confi_DORIS_USO_correction', confi_DORIS_USO_correction,...
                            'confi_complex_cal1_correction_from_IPF_DB', ins_star_tracker_2,...
                            'confi_TRK_echo_error', ins_star_tracker_3,...
                            'confi_echo_rx1_error',ins_reserved4,...                            
                            'confi_echo_rx2_error', confi_echo_rx2_error,...
                            'confi_NMP_inconsistency', confi_NMP_inconsistency,...
                            'confi_azimuth_cal_missing', confi_azimuth_cal_missing,...
                            'confi_azimuth_cal_from_IPF_DB', confi_azimuth_cal_from_IPF_DB,...
                            'confi_range_window_cal_function_missing',confi_range_window_cal_function_missing,...
                            'confi_range_window_cal_function_from_IPF_DB', confi_range_window_cal_function_from_IPF_DB,...
                            'confi_reserved', confi_reserved,...
                            'confi_cal2_correction_missing', confi_cal2_correction_missing,...
                            'confi_cal2_correction_from_IPF_DB', confi_cal2_correction_from_IPF_DB,...
                            'confi_power_scaling_error_LRM_only', confi_power_scaling_error_LRM_only,...
                            'confi_attitude_correction_missing', confi_attitude_correction_missing,...
                            'confi_attitude_interpolation_error',confi_phase_perturbation_correction,...
                            'confi_reserved2',confi_reserved2);

        time_orbit_group(j)=struct(    'data_record_time', data_record_time,...
                                'USO_correction', USO_correction,...
                                'modeID', modeID,...
                                'ssc', ssc,...
                                'instrument_configuration', instrument_configuration,...
                                'burst_counter',burst_counter,...
                                'latitude',latitude,...
                                'longitude',longitude,...
                                'altitude_COG',altitude_COG,...
                                'instantaneous_altitude_rate',instantaneous_altitude_rate,...
                                'satellite_velocity',satellite_velocity,...
                                'real_beam_direction',real_beam_direction,...
                                'inferometer_baseline',inferometer_baseline,...
                                'str_id'  ,str_id,...
                                'bench_roll',bench_roll,...
                                'bench_pitch' ,bench_pitch,...
                                'bench_yaw',bench_yaw,...
                                'L1B_measurement_confidence_flag',L1B_measurement_confidence_flag);

    end
        
               
    
   
    
    %% read Measurements Group %
    
    
    for j=1:20
        
                
        %15 Window Delay (2way) uncorrected for instrument delays 10-12 s 8 sll
        window_delay = fread(fid,1,'int64');
        %16 Ho Initial Height Word 48.8 ps 4 sl (see note 2)
        ho_initial_height  = fread(fid,1,'int32');
        %17 HPR Height Rate 3.05 ps/rc 4 sl (see note 2)
        HPR_hight_rate = fread(fid,1,'int32');
        %18 LAI 12.5 ns 4 sl (see note 2)
        LAI = fread(fid,1,'int32');
        %19 FAI 12.5/256 ns 4 sl (see note 2)
        FAI = fread(fid,1,'int32');
        %20 AGC_1 (not corrected) dB/100 4 sl (see note 3)
        AGC_1 = fread(fid,1,'int32');
        %21 AGC_2 (not corrected) dB/100 4 sl (see note 3)
        AGC_2 = fread(fid,1,'int32');
        %22 Total Fixed Gain Rx 1 dB/100 4 sl (see note 3)
        tot_fixed_gain_1 = fread(fid,1,'int32');
        %23 Total Fixed Gain Rx 2 dB/100 4 sl (see note 3)
        tot_fixed_gain_2 = fread(fid,1,'int32');
        %24 Transmit Power Micro-Watts 4 sl
        transmit_power = fread(fid,1,'int32');
        %25 Doppler range correction (Radial component)mm 4 sl
        doppler_range_correction = fread(fid,1,'int32');
        %26 Instrument Range Correction tx-rx antenna mm 4 sl
        instrument_range_correction_tx_rx = fread(fid,1,'int32');
        %27 Instrument Range Correction mm 4 sl
        instrument_range_correction = fread(fid,1,'int32');
        %28 Instrument Sigma 0 correction, tx-rx antenna dB/100 4 sl (see note
        instrument_sigma0_correction_tx_rx = fread(fid,1,'int32');
        %29 Instrument Sigma 0 correction rx only antenna dB/100 4 sl (see note
        instrument_sigma0_correction_rx = fread(fid,1,'int32');
        %30 Internal Phase Correction Microradians 4 sl
        internal_phase_correction = fread(fid,1,'int32');
        %31 External Phase Correction Microradians 4 sl
        external_phase_correction = fread(fid,1,'int32');
        %32 Noise power measurement dB/100 4 sl (see note
        noise_power = fread(fid,1,'int32');
        %33 Phase Slope Correction Microradians 4 sl (see note
        phase_slope_correction = fread(fid,1,'int32');
        %34 spares 4*1 uc
        fread(fid,4,'uchar');
        measurements_group(j)=struct(   'window_delay',window_delay,...
                                        'ho_initial_height',ho_initial_height,...
                                        'HPR_hight_rate',HPR_hight_rate,...
                                        'LAI',LAI,...
                                        'FAI',FAI,...
                                        'AGC_1',AGC_1,...
                                        'AGC_2',AGC_2,...
                                        'tot_fixed_gain_1',tot_fixed_gain_1,...
                                        'tot_fixed_gain_2',tot_fixed_gain_2,...
                                        'transmit_power',transmit_power,...
                                        'doppler_range_correction',doppler_range_correction,...
                                        'instrument_range_correction_tx_rx',instrument_range_correction_tx_rx,...
                                        'instrument_range_correction',instrument_range_correction,...
                                        'instrument_sigma0_correction_tx_rx',instrument_sigma0_correction_tx_rx,...
                                        'instrument_sigma0_correction_rx',instrument_sigma0_correction_rx,...
                                        'internal_phase_correction',internal_phase_correction,...
                                        'external_phase_correction',external_phase_correction,...
                                        'phase_slope_correction',phase_slope_correction,...
                                        'noise_power',noise_power);
                                        
           
    end
    
    %%  read Corrections Group %
        
    %35 Dry Tropospheric Correction mm 4 sl
    dry_tropo_correction = fread(fid,1,'int32');
    %36 Wet Tropospheric Correction mm 4 sl
    wet_tropo_correction = fread(fid,1,'int32');
    %37 Inverse Barometric Correction mm 4 sl
    inverse_baro_correction = fread(fid,1,'int32');
    %38 Dynamic Atmospheric Correction mm 4 sl
    Dynamic_atmospheric_correction = fread(fid,1,'int32');        
    %39 GIM Ionospheric Correction mm 4 sl
    GIM_iono_correction = fread(fid,1,'int32');
    %40 Model Ionospheric Correction mm 4 sl
    model_iono_correction = fread(fid,1,'int32');
    %41 Ocean Equilibrium Tide mm 4 sl
    ocean_equilibrium_tide = fread(fid,1,'int32');
    %42 Long Period Tide Height 4 sl
    long_period_tide_height = fread(fid,1,'int32'); 
    %43 Ocean Loading Tide mm 4 sl
    ocean_loading_tide = fread(fid,1,'int32');
    %44 Solid Earth Tide mm 4 sl
    solid_earth_tide = fread(fid,1,'int32');
    %45 Geocentric Polar Tide mm 4 sl
    geocentric_polar_tide = fread(fid,1,'int32');
    %46 Surface type flag 4 ul
    surface_type_flag = fread(fid,1,'uint32');
    %47 spares 4*1 uc
    fread(fid,4,'uchar');
    %48 Correction status flags 4 ul (see table 2.3.3-5)
    %correction_status_flags = fread(fid,1,'uint32');

    correction_status_flags(1) = fread(fid,1,'ubit1');
    correction_status_flags(2) = fread(fid,1,'ubit1');
    correction_status_flags(3) = fread(fid,1,'ubit1');
    correction_status_flags(4) = fread(fid,1,'ubit1');
    correction_status_flags(5) = fread(fid,1,'ubit1');
    correction_status_flags(6) = fread(fid,1,'ubit1');
    correction_status_flags(7) = fread(fid,1,'ubit1');
    correction_status_flags(8) = fread(fid,1,'ubit1');
    correction_status_flags(9) = fread(fid,1,'ubit1');
    correction_status_flags(10) = fread(fid,1,'ubit1');
    correction_status_flags(11) = fread(fid,1,'ubit1');
    correction_status_flags(12) = fread(fid,1,'ubit1');
    correction_status_flags(13) = fread(fid,1,'ubit20');


    %49 Correction error flags 4 ul (see table 2.3.3-6)
    %correction_error_flags = fread(fid,1,'uint32');

    correction_error_flags(1) = fread(fid,1,'ubit1');
    correction_error_flags(2) = fread(fid,1,'ubit1');
    correction_error_flags(3) = fread(fid,1,'ubit1');
    correction_error_flags(4) = fread(fid,1,'ubit1');
    correction_error_flags(5) = fread(fid,1,'ubit1');
    correction_error_flags(6) = fread(fid,1,'ubit1');
    correction_error_flags(7) = fread(fid,1,'ubit1');
    correction_error_flags(8) = fread(fid,1,'ubit1');
    correction_error_flags(9) = fread(fid,1,'ubit1');
    correction_error_flags(10) = fread(fid,1,'ubit1');
    correction_error_flags(11) = fread(fid,1,'ubit1');
    correction_error_flags(12) = fread(fid,1,'ubit1');
    correction_error_flags(13) = fread(fid,1,'ubit20');


    %50 Spare 4*1 uc
    fread(fid,4,'uchar');


    corrections_group = struct( 'dry_tropo_correction',dry_tropo_correction,...
                                'wet_tropo_correction',wet_tropo_correction,...
                                'inverse_baro_correction',inverse_baro_correction,...
                                'Dynamic_atmospheric_correction',Dynamic_atmospheric_correction,...
                                'GIM_iono_correction',GIM_iono_correction,...
                                'model_iono_correction',model_iono_correction,...
                                'ocean_equilibrium_tide',ocean_equilibrium_tide,...
                                'long_period_tide_height',long_period_tide_height,...
                                'ocean_loading_tide',ocean_loading_tide,...
                                'solid_earth_tide',solid_earth_tide,...
                                'geocentric_polar_tide',geocentric_polar_tide,...
                                'surface_type_flag',surface_type_flag,...
                                'correction_status_flags',correction_status_flags,...
                                'correction_error_flags',correction_error_flags);

    %   disp(['Corrections Group ' num2str(i)]);                    
    %   disp(['dry_tropo_correction ' num2str(dry_tropo_correction)]);
    % 	disp(['wet_tropo_correction ' num2str(wet_tropo_correction)]);
    % 	disp(['inverse_baro_correction ' num2str(inverse_baro_correction)]);
    % 	disp(['DORIS_iono_correction ' num2str(DORIS_iono_correction)]);
    % 	disp(['model_iono_correction ' num2str(model_iono_correction)]);
    % 	disp(['ocean_tide_1_mm ' num2str(ocean_tide_1_mm)]);
    % 	disp(['ocean_tide_2_mm ' num2str(ocean_tide_2_mm)]);
    % 	disp(['ocean_loading_tide ' num2str(ocean_loading_tide)]);
    % 	disp(['solid_earth_tide ' num2str(solid_earth_tide)]);
    % 	disp(['geocentric_polar_tide ' num2str(geocentric_polar_tide)]);
    % 	disp(['surface_type_flag ' num2str(surface_type_flag)]);
    % 	disp(['correction_status_flags ' num2str(correction_status_flags)]);
    % 	disp(['correction_error_flags ' num2str(correction_error_flags)]);


     %% read Average Waveforms Group (LRM) %
   
   if(strcmp(siralMode,'LRM'))
            % 51 Data Record Time (MDSR Time Stamp) TAI 12 sl+2*ul
            day = fread(fid,1,'int32');
            second = fread(fid,1,'uint32');
            microsecond = fread(fid, 1, 'uint32');
            data_record_time = struct(  'day', day,...
                                        'second', second, ...
                                        'microsecond',microsecond);

            % 52 Latitude of measurement 10-1 ìdeg 4 sl (see note 1)
            latitude = fread(fid,1,'int32');

            % 53 Longitude of measurement 10-1 ìdeg 4 sl (see note 1)
            longitude = fread(fid,1,'int32');

            % 54 Altitude of COG above reference ellipsoid(interpolated value) mm 4 sl
            altitude = fread(fid,1,'int32');

            % 55 Window Delay (2way) corrected for instrument delays 10-12 s 8 sll
            window_delay = fread(fid,1,'int64');

            % 56 1 Hz Averaged Power Echo Waveform[128] Scaled 512*2 us
            n_samples=128;
            for iSample=1:n_samples
               averaged_echo_waveform(iSample) = fread(fid,1,'uint16');       
            end

            % 57 Echo Scale Factor (to scale echo to watts) - 4 sl
            echo_scale_factor = fread(fid,1,'int32');

            % 58 Echo Scale Power (a power of 2) 4 sl
            echo_scale_power = fread(fid,1,'int32');

            % 59 Number of echoes averaged - 2 us
            number_echoes = fread(fid,1,'uint16');
            % 60 Flags - 2 us
            flags = fread(fid,1,'uint16');
            average_waveform_group = struct(    'data_record_time',data_record_time,...
                                                'latitude',latitude,...
                                                'longitude',longitude,...
                                                'altitude',altitude,...
                                            'window_delay',window_delay,...
                                            'averaged_echo_waveform',averaged_echo_waveform,...
                                            'echo_scale_factor',echo_scale_factor,...
                                            'echo_scale_power',echo_scale_power,...
                                            'number_echoes',number_echoes,...
                                            'flags',flags);
        %	disp(['Average Waveform ' num2str(i)]);                                    
        % 	disp(['data_record_time ' num2str(data_record_time)]);
        % 	disp(['latitude ' num2str(latitude)]);
        % 	disp(['longitude ' num2str(longitude)]);
        % 	disp(['altitude ' num2str(altitude)]);
        % 	disp(['window_delay ' num2str(window_delay)]);
        % 	disp(['averaged_echo_waveform ' num2str(averaged_echo_waveform)]);
        % 	disp(['echo_scale_factor ' num2str(echo_scale_factor)]);
        % 	disp(['echo_scale_power ' num2str(echo_scale_power)]);
        % 	disp(['number_echoes ' num2str(number_echoes)]);
        % 	disp(['flags ' num2str(flags)]);
   end
    
   %% read Average Waveforms Group (SAR) %
   
   if(strcmp(siralMode,'SAR'))
            % 51 Data Record Time (MDSR Time Stamp) TAI 12 sl+2*ul
            day = fread(fid,1,'int32');
            second = fread(fid,1,'uint32');
            microsecond = fread(fid, 1, 'uint32');
            data_record_time = struct(  'day', day,...
                                        'second', second, ...
                                        'microsecond',microsecond);

            % 52 Latitude of measurement 10-1 ìdeg 4 sl (see note 1)
            latitude = fread(fid,1,'int32');

            % 53 Longitude of measurement 10-1 ìdeg 4 sl (see note 1)
            longitude = fread(fid,1,'int32');

            % 54 Altitude of COG above reference ellipsoid(interpolated value) mm 4 sl
            altitude = fread(fid,1,'int32');

            % 55 Window Delay (2way) corrected for instrument delays 10-12 s 8 sll
            window_delay = fread(fid,1,'int64');

            % 56 1 Hz Averaged Power Echo Waveform[128] Scaled 512*2 us
            n_samples=128;
            for iSample=1:n_samples
               averaged_echo_waveform(iSample) = fread(fid,1,'uint16');       
            end

            % 57 Echo Scale Factor (to scale echo to watts) - 4 sl
            echo_scale_factor = fread(fid,1,'int32');

            % 58 Echo Scale Power (a power of 2) 4 sl
            echo_scale_power = fread(fid,1,'int32');

            % 59 Number of echoes averaged - 2 us
            number_echoes = fread(fid,1,'uint16');
            % 60 Flags - 2 us
            flags = fread(fid,1,'uint16');
            average_waveform_group = struct(    'data_record_time',data_record_time,...
                                                'latitude',latitude,...
                                                'longitude',longitude,...
                                                'altitude',altitude,...
                                            'window_delay',window_delay,...
                                            'averaged_echo_waveform',averaged_echo_waveform,...
                                            'echo_scale_factor',echo_scale_factor,...
                                            'echo_scale_power',echo_scale_power,...
                                            'number_echoes',number_echoes,...
                                            'flags',flags);
        %	disp(['Average Waveform ' num2str(i)]);                                    
        % 	disp(['data_record_time ' num2str(data_record_time)]);
        % 	disp(['latitude ' num2str(latitude)]);
        % 	disp(['longitude ' num2str(longitude)]);
        % 	disp(['altitude ' num2str(altitude)]);
        % 	disp(['window_delay ' num2str(window_delay)]);
        % 	disp(['averaged_echo_waveform ' num2str(averaged_echo_waveform)]);
        % 	disp(['echo_scale_factor ' num2str(echo_scale_factor)]);
        % 	disp(['echo_scale_power ' num2str(echo_scale_power)]);
        % 	disp(['number_echoes ' num2str(number_echoes)]);
        % 	disp(['flags ' num2str(flags)]);
   end
        
   %% read Average Waveforms Group (SARin) %
   
   if(strcmp(siralMode,'SARin'))
            % 61 Data Record Time (MDSR Time Stamp) TAI 12 sl+2*ul
            day = fread(fid,1,'int32');
            second = fread(fid,1,'uint32');
            microsecond = fread(fid, 1, 'uint32');
            data_record_time = struct(  'day', day,...
                                        'second', second, ...
                                        'microsecond',microsecond);

            % 62 Latitude of measurement 10-1 ìdeg 4 sl (see note 1)
            latitude = fread(fid,1,'int32');

            % 63 Longitude of measurement 10-1 ìdeg 4 sl (see note 1)
            longitude = fread(fid,1,'int32');

            % 64 Altitude of COG above reference ellipsoid(interpolated value) mm 4 sl
            altitude = fread(fid,1,'int32');

            % 65 Window Delay (2way) corrected for instrument delays 10-12 s 8 sll
            window_delay = fread(fid,1,'int64');

            % 66 1 Hz Averaged Power Echo Waveform[512] Scaled 512*2 us
            n_samples=512;
            for iSample=1:n_samples
               averaged_echo_waveform(iSample) = fread(fid,1,'uint16');       
            end

            % 67 Echo Scale Factor (to scale echo to watts) - 4 sl
            echo_scale_factor = fread(fid,1,'int32');

            % 68 Echo Scale Power (a power of 2) 4 sl
            echo_scale_power = fread(fid,1,'int32');

            % 69 Number of echoes averaged - 2 us
            number_echoes = fread(fid,1,'uint16');
            % 70 Flags - 2 us
            flags = fread(fid,1,'uint16');
            average_waveform_group = struct(    'data_record_time',data_record_time,...
                                                'latitude',latitude,...
                                                'longitude',longitude,...
                                                'altitude',altitude,...
                                            'window_delay',window_delay,...
                                            'averaged_echo_waveform',averaged_echo_waveform,...
                                            'echo_scale_factor',echo_scale_factor,...
                                            'echo_scale_power',echo_scale_power,...
                                            'number_echoes',number_echoes,...
                                            'flags',flags);
        %	disp(['Average Waveform ' num2str(i)]);                                    
        % 	disp(['data_record_time ' num2str(data_record_time)]);
        % 	disp(['latitude ' num2str(latitude)]);
        % 	disp(['longitude ' num2str(longitude)]);
        % 	disp(['altitude ' num2str(altitude)]);
        % 	disp(['window_delay ' num2str(window_delay)]);
        % 	disp(['averaged_echo_waveform ' num2str(averaged_echo_waveform)]);
        % 	disp(['echo_scale_factor ' num2str(echo_scale_factor)]);
        % 	disp(['echo_scale_power ' num2str(echo_scale_power)]);
        % 	disp(['number_echoes ' num2str(number_echoes)]);
        % 	disp(['flags ' num2str(flags)]);
   end
   
   

    
     %% read Waveform Group (LRM) %

    if(strcmp(siralMode,'LRM'))
        
        for j=1:20
            % 71 Averaged Power Echo Waveform[128] Scaled 128*2 us

%             if (record_num==46 && j==20)
%                 j=j;
%             end

            for iSample=1:128
                averaged_power_echo_waveform(iSample) = fread(fid,1,'uint16');
            end

            % 72 Echo Scale Factor (to scale echo to watts) - 4 sl
            echo_scale_factor = fread(fid,1,'int32');
            % 73 Echo Scale Power (a power of 2) 4 sl
            echo_scale_power = fread(fid,1,'int32');
            % 74 Number of echoes averaged - 2 us
            number_echoes = fread(fid,1,'uint16');
            % 75 Flags - 2 us
            flags = fread(fid,1,'uint16');
            
            waveform_group(j) = struct( 'averaged_power_echo_waveform',averaged_power_echo_waveform,...
                                        'echo_scale_factor',echo_scale_factor,...
                                        'echo_scale_power',echo_scale_power,...
                                        'number_echoes',number_echoes,...
                                        'flags',flags);
                                        
    %       disp(['Waveform ' num2str(j) '/' num2str(i)]);                                
    % 		disp(['averaged_power_echo_waveform ' num2str(averaged_power_echo_waveform)]);
    % 		disp(['echo_scale_factor ' num2str(echo_scale_factor)]);
    % 		disp(['echo_scale_power ' num2str(echo_scale_power)]);
    % 		disp(['number_echoes ' num2str(number_echoes)]);
    % 		disp(['flags ' num2str(flags)]);
    % 		disp(['beam_behavour_parameter ' num2str(beam_behavour_parameter)]);
    % 		disp(['coherence ' num2str(coherence)]);
    % 		disp(['phase_difference ' num2str(phase_difference)]);
            
        end
    end
     %% read Waveform Group (SAR) %
    
    if(strcmp(siralMode,'SAR'))
        
        for j=1:20
            % 76 Averaged Power Echo Waveform[128] Scaled 128*2 us

%             if (record_num==46 && j==20)
%                 j=j;
%             end

           for iSample=1:256
                averaged_power_echo_waveform(iSample) = fread(fid,1,'uint16');
            end

            % 77 Echo Scale Factor (to scale echo to watts) - 4 sl
            echo_scale_factor = fread(fid,1,'int32');
            % 78 Echo Scale Power (a power of 2) 4 sl
            echo_scale_power = fread(fid,1,'int32');
            % 79 Number of echoes averaged - 2 us
            number_echoes = fread(fid,1,'uint16');
            % 80 Flags - 2 us
            flags = fread(fid,1,'uint16');
            % 81 Beam behaviour parameter [50] 50*2 ss (see table 2.3.5-4)
            for iBeam=1:50
                beam_behavour_parameter(iBeam)=fread(fid,1,'int16');
            end
            waveform_group(j) = struct( 'averaged_power_echo_waveform',averaged_power_echo_waveform,...
                                        'echo_scale_factor',echo_scale_factor,...
                                        'echo_scale_power',echo_scale_power,...
                                        'number_echoes',number_echoes,...
                                        'flags',flags,...
                                        'beam_behavour_parameter',beam_behavour_parameter);
    %       disp(['Waveform ' num2str(j) '/' num2str(i)]);                                
    % 		disp(['averaged_power_echo_waveform ' num2str(averaged_power_echo_waveform)]);
    % 		disp(['echo_scale_factor ' num2str(echo_scale_factor)]);
    % 		disp(['echo_scale_power ' num2str(echo_scale_power)]);
    % 		disp(['number_echoes ' num2str(number_echoes)]);
    % 		disp(['flags ' num2str(flags)]);
    % 		disp(['beam_behavour_parameter ' num2str(beam_behavour_parameter)]);
    % 		disp(['coherence ' num2str(coherence)]);
    % 		disp(['phase_difference ' num2str(phase_difference)]);
            
        end
     end
%     disp(iRecord);
        
    
    %% read Waveform Group (SARin) %
    
    if(strcmp(siralMode,'SARin'))
       
        for j=1:20
            % 82 Averaged Power Echo Waveform[512] Scaled 512*2 us

%             if (record_num==46 && j==20)
%                 j=j;
%             end

            for iSample=1:1024
                averaged_power_echo_waveform(iSample) = fread(fid,1,'uint16');
            end

            % 83 Echo Scale Factor (to scale echo to watts) - 4 sl
            echo_scale_factor = fread(fid,1,'int32');
            % 84 Echo Scale Power (a power of 2) 4 sl
            echo_scale_power = fread(fid,1,'int32');
            % 85 Number of echoes averaged - 2 us
            number_echoes = fread(fid,1,'uint16');
            % 86 Flags - 2 us
            flags = fread(fid,1,'uint16');
            % 87 Beam behaviour parameter [50] 50*2 ss (see table 2.3.5-4)
            for iBeam=1:50
                beam_behavour_parameter(iBeam)=fread(fid,1,'int16');
            end
            % 88 Coherence [512] 1/1000 512*2 us
            for iSample=1:1024
                coherence(iSample) = fread(fid,1,'uint16');
            end
            % 89 Phase difference [512] microrad 512*4 sl
            for iSample=1:1024
                phase_difference(iSample) = fread(fid,1,'int32');
            end
            waveform_group(j) = struct( 'averaged_power_echo_waveform',averaged_power_echo_waveform,...
                                        'echo_scale_factor',echo_scale_factor,...
                                        'echo_scale_power',echo_scale_power,...
                                        'number_echoes',number_echoes,...
                                        'flags',flags,...
                                        'beam_behavour_parameter',beam_behavour_parameter,...
                                        'coherence',coherence,...
                                        'phase_difference',phase_difference);
    %       disp(['Waveform ' num2str(j) '/' num2str(i)]);                                
    % 		disp(['averaged_power_echo_waveform ' num2str(averaged_power_echo_waveform)]);
    % 		disp(['echo_scale_factor ' num2str(echo_scale_factor)]);
    % 		disp(['echo_scale_power ' num2str(echo_scale_power)]);
    % 		disp(['number_echoes ' num2str(number_echoes)]);
    % 		disp(['flags ' num2str(flags)]);
    % 		disp(['beam_behavour_parameter ' num2str(beam_behavour_parameter)]);
    % 		disp(['coherence ' num2str(coherence)]);
    % 		disp(['phase_difference ' num2str(phase_difference)]);
            
        end
    end

    %%  OUTPUT  --%
    %--------------%
    
        l1b_dsr = struct(   'time_orbit_group',time_orbit_group,...
                            'measurements_group',measurements_group,...
                            'corrections_group',corrections_group,...
                            'average_waveform_group',average_waveform_group,...
                            'waveform_group',waveform_group);

        clear time_orbit_group
        clear measurements_group
        clear corrections_group
        clear average_waveform_group
        clear waveform_group

        l1b_ds(iRecord) = l1b_dsr;
        out_L1B=l1b_ds;
    
end
fclose(fid);
% averaged_echo_waveform=zeros(360,128);
% for kkk=1:record_num
%     averaged_echo_waveform(kkk,:)=out_L1b(1,kkk).average_waveform_group.averaged_echo_waveform(:);
%     
%     
% end
% meshc(averaged_echo_waveform);hold on

    


%% Plotting

    if(KML)
       bandwidth_value = 320*1e6;
       c=2.99792458e8;
       Waveform  = zeros(length(out_L1B),n_samples);
       Time      = zeros(1,length(out_L1B));
       Latitude  = zeros(1,length(out_L1B));
       Longitude = zeros(1,length(out_L1B));
       altitude = zeros(1,length(out_L1B));
       Latitude_DB  = zeros(1,length(out_L1B)*20);
       Longitude_DB = zeros(1,length(out_L1B)*20);
       altitude_DB = zeros(1,length(out_L1B)*20);
       window_delay = zeros(1,length(out_L1B));
       max_sample = zeros(1,length(out_L1B));
       window_delay_DB = zeros(1,length(out_L1B)*20);
       max_sample_DB = zeros(1,length(out_L1B)*20);
       z=zeros(1,length(out_L1B));
       bandwidth_value = 320*1e6;
       c=2.99792458e8;
       for iRecord=1:length(out_L1B)
%            customText = 'Loop...';
%             percentageDone =  iRecord / length(out_L1B);
%             stopBar= progressbar(percentageDone, 0, customText);
%             if (stopBar) 
%                 break; 
%             end
%            Waveform(iRecord,:)=(out_L1B(iRecord).average_waveform_group.averaged_echo_waveform);
%            figure; plot(Waveform(iRecord,:));
%            set(gca,'FontSize',20) % --> eixos
%            xlab = get(gca,'XLabel'); set(xlab,'String','Samples','FontSize',18,'FontName','Arial')
%            ylab = get(gca,'YLabel'); set(ylab,'String','Power','FontSize',18,'FontName','Arial')
%            saveas(gcf, sprintf('Waveform%d.png',iRecord) , 'png');
%            close all
%            fclose all 
% %            Time(iRecord)=946684800+(out_L1B(iRecord).average_waveform_group.data_record_time.day)*3600*24 ... %946684800 UTC second of starting point 01/01/2000 0:00:00 
% %                +(out_L1B(iRecord).average_waveform_group.data_record_time.second)+ ... 
% %                1e-6*(out_L1B(iRecord).average_waveform_group.data_record_time.microsecond);
%            Latitude(iRecord)=(out_L1B(iRecord).average_waveform_group.latitude*1e-7);
%            Longitude(iRecord)=(out_L1B(iRecord).average_waveform_group.longitude*1e-7);
%            altitude(iRecord)=out_L1B(iRecord).average_waveform_group.altitude;
%            window_delay(iRecord)=out_L1B(iRecord).average_waveform_group.window_delay;
           
           dry(iRecord) = out_L1B(iRecord).corrections_group.dry_tropo_correction;
%            [~,max_sample(iRecord)]=max(Waveform(iRecord,:));
           
           for iDB=1:20
               
            
               window_delay_DB((iRecord-1)*20+iDB)=out_L1B(iRecord).measurements_group(iDB).window_delay;
%                doppler_range_correction((iRecord-1)*20+iDB)=out_L1B(iRecord).measurements_group(iDB).doppler_range_correction;
               Waveform_DB((iRecord-1)*20+iDB,:)=(out_L1B(iRecord).waveform_group(iDB).averaged_power_echo_waveform);
               Latitude_DB((iRecord-1)*20+iDB)=(out_L1B(iRecord).time_orbit_group(iDB).latitude*1e-7);
               Longitude_DB((iRecord-1)*20+iDB)=(out_L1B(iRecord).time_orbit_group(iDB).longitude*1e-7);
               altitude_DB((iRecord-1)*20+iDB)=out_L1B(iRecord).time_orbit_group(iDB).altitude_COG;
               [~,max_sample_DB((iRecord-1)*20+iDB)]=max(Waveform_DB((iRecord-1)*20+iDB,:));
               fine_height((iRecord-1)*20+iDB) = (altitude_DB((iRecord-1)*20+iDB)-(window_delay_DB((iRecord-1)*20+iDB)*1e-9 +(max_sample_DB((iRecord-1)*20+iDB)-((n_samples)/2))*(1/(bandwidth_value)))*c/2) * 1e-3; %height in meters 
                Phase((iRecord-1)*20+iDB,:)=(out_L1B(iRecord).waveform_group(iDB).phase_difference);
                Coherence((iRecord-1)*20+iDB,:)=(out_L1B(iRecord).waveform_group(iDB).coherence);
%                inferometer_baseline((iRecord-1)*20+iDB,:) = out_L1B(iRecord).time_orbit_group(iDB).inferometer_baseline;
%                fprintf(fidTXT,'%d  %d  %d \n',Latitude_DB,Longitude_DB,fine_height);
%                roll((iRecord-1)*20+iDB)  = atan(inferometer_baseline((iRecord-1)*20+iDB,1)/inferometer_baseline((iRecord-1)*20+iDB,3)); % [rad]
               
           end
       end
%        fprintf(fidTXT,';\n');
       
%     for k=1:floor(size(Waveform_DB,1)/20)
% 
%         Waveform(k,:)=mean(Waveform_DB((k-1)*20+1:k*20,:),1);
%         window_delay(k)=mean(window_delay_DB((k-1)*20+1:k*20));
%         altitude(k)=mean(altitude_DB((k-1)*20+1:k*20));
%         Latitude(k)=mean(Latitude_DB((k-1)*20+1:k*20));
%         Longitude(k)=mean(Longitude_DB((k-1)*20+1:k*20));
%         
%         [max_val_L1_1Hz(k),max_pos_CR_L1_1Hz(k)]=max(Waveform(k,:));    
%         [retrack_val_L1_1Hz(k),retrack_pos_CR_L1_1Hz(k)]=min(abs(Waveform(k,1:max_pos_CR_L1_1Hz(k))-0.87*max_val_L1_1Hz(k))); % L2 retracker point (Cristina)
%         coarse_height(k)=(altitude(k)*1e-3-(window_delay(k)*1e-12+(retrack_pos_CR_L1_1Hz(k)-(n_samples)/2-1)/bandwidth_value)*c/2);
%     end
%        n_samples
%        plot(doppler_range_correction)
       
%        coarse_height = (altitude-(window_delay*1e-9 +(max_sample-((n_samples)/2))*(1/(bandwidth_value)))*c/2) * 1e-3; %height in meters 
       
        lla2kml(filename,Latitude_DB,Longitude_DB,fine_height,'.');
 
 

%               
%        latAxis = [Latitude(1):(Latitude(end)-Latitude(1))/(length(out_L1B)*20-1):Latitude(end)];
%        figure;plot((-roll+asin(Phase(1220,:)*1e-6*22.084/(2*pi*1167))*180/pi),'.r');
%        figure;meshc(1:n_samples,latAxis,Waveform_DB);
% %        AoA = asin(lambda*phase/(2*pi*B)
% %        meshc(1:512,latAxis(117*20:119*20),asin (22.0840 * (Phase(118*20:119*20,:)./(2 * pi * 1000000 * 1.1676e+003))*180/pi));
% 
%        set(get(gca,'XLabel'),'String','Samples');
%        set(get(gca,'YLabel'),'String','Latitude [deg]');
%        set(get(gca,'ZLabel'),'String','');
%        view((Latitude(1)-Latitude(end)),85);
%        saveas(gcf, sprintf('%s_waveforms.fig',filename) , 'fig')
%        saveas(gcf, sprintf('%s_waveforms.jpg',filename) , 'jpg')
%        
%        fclose all
%        
%         figure
%         usamap([1e-6*Latitude(end) 1e-6*Latitude(1)], [1e-6*Longitude(end) 1e-6*Longitude(1)])
%         geoshow(1e-6*Latitude, 1e-6*Longitude, z, 'DisplayType', 'surface')
%         demcmap(z)
%         daspectm('m',1)
%         view(3)
    end
end



