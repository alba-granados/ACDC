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
% READFBR: function that reads the FBR data set from input filename,
% assuming SARIN data records
%
% Calling
%   fbr_ds = readFBR( filename, headers )
%
% Inputs
%   filename: input SARIN FBR file
%   headers:   if false --> without header
%              if true --> with header
%
% Output
%   fbr_ds           : data contained in the file 
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
%           Daniel Martínez / Pildo Labs
%           Josep Montolio / Pildo Labs
%           Pablo García / isardSAT
%           Albert Garcia /isardSAT
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Mònica Roca / isardSAT (26/05/11)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
% Now available for SAR and SARin
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [fbr_ds,record_num] = auto_readFBR(filename, headers,mode)

end_of_file = 0;
record_num=0;

if ( headers )
%     disp(['reading FBR from ' filename ' WITH headers']);
    sph = readSPH(filename);
    % FBR data set must be the first one
    fbr_dsd = sph.dsds(1);
    if(~strfind(fbr_dsd.ds_name,'SIR_FBR'))
        disp('no SIR_FBR_SARIN data set found');
        return;
    else
%         disp('found SIR_FBR_SARIN data set');
%         disp(fbr_dsd);
    end

    ds_offset = fbr_dsd.ds_offset;
    ds_size = fbr_dsd.ds_size;
    num_dsr = fbr_dsd.num_dsr;
    dsr_size = fbr_dsd.dsr_size;
else
    disp(['reading FBR from ' filename ' WITHOUT headers'] );
    
    ds_offset = 0;        
end

% read the dsrs
fid = fopen(filename,'r','b');
fseek(fid, ds_offset, 0);

%iterate available dsr's

for record_num=1:num_dsr
% while (end_of_file==0)
        
%     record_num = record_num+1; %One record is added to the record counter
    
    %--------------------------------%
    %--  read Time and Orbit Group --%
    %--------------------------------%
    
    
    for j=1:20
		%1 Data Record Time (MDSR TimeStamp) TAI 12 sl+2*ul
        day = fread(fid,1,'int32');
        
        if(isempty(day)==1)%If we don't find anything inside means that the file has end
            end_of_file=1;
%             record_num = record_num-1;
        else
        
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
            %13 Interferometer baseline vector[3](in CRF)micrometers 3*4 sl
            inferometer_baseline(1) = fread(fid,1,'int32');
            inferometer_baseline(2) = fread(fid,1,'int32');
            inferometer_baseline(3) = fread(fid,1,'int32');
            %14 FBR Measurement ConfidenceData (flag word)4 u
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
                confi_attitude_interpolation_error = fread(fid,1,'ubit1');
                confi_reserved2 = fread(fid,1,'ubit1');
                confi_phase_perturbation = fread(fid,1,'ubit1');

                FBR_measurement_confidence_flag =struct(    'confi_block_degraded', confi_block_degraded,...
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
                                    'confi_attitude_interpolation_error',confi_attitude_interpolation_error,...
                                    'confi_reserved2',confi_reserved2,...
                                    'confi_phase_perturbation',confi_phase_perturbation);

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
                                        'FBR_measurement_confidence_flag',FBR_measurement_confidence_flag);



        %       disp(['Time and Orbit Group ' num2str(j) '/' num2str(i)]);
        %       disp(['data_record_time' num2str(data_record_time)]);
        % 		disp(['USO_correction ' num2str(USO_correction)]);
        % 		disp(['modeID ' num2str(modeID)]);
        % 		disp(['ssc ' num2str(ssc)]);
        % 		disp(['instrument_configuration ' num2str(instrument_configuration)]);
        % 		disp(['burst_counter ' num2str(burst_counter)]);
        % 		disp(['latitude ' num2str(latitude)]);
        % 		disp(['longitude ' num2str(longitude)]);
        % 		disp(['altitude_COG ' num2str(altitude_COG)]);
        % 		disp(['instantaneous_altitude_rate ' num2str(instantaneous_altitude_rate)]);
        % 		disp(['satellite_velocity ' num2str(satellite_velocity)]);
        % 		disp(['real_beam_direction ' num2str(real_beam_direction)]);
        % 		disp(['inferometer_baseline ' num2str(inferometer_baseline)]);
        % 		disp(['FBR_measurement_confidence_flag ' num2str(FBR_measurement_confidence_flag)]);       
        end
%         j = j+1;
    end
    
    %-----------------------------%
    %-- read Measurements Group --%
    %-----------------------------%
          
    
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
                                        'noise_power',noise_power);
                                            
%       disp(['Measurements Group ' num2str(j) '/' num2str(i)]);
% 		disp(['window_delay ' num2str(window_delay)]);
% 		disp(['ho_initial_height ' num2str(ho_initial_height)]);
% 		disp(['HPR_hight_rate ' num2str(HPR_hight_rate)]);
% 		disp(['LAI ' num2str(LAI)]);
% 		disp(['FAI ' num2str(FAI)]);
% 		disp(['AGC_1 ' num2str(AGC_1)]);
% 		disp(['AGC_2 ' num2str(AGC_2)]);
% 		disp(['tot_fixed_gain_1 ' num2str(tot_fixed_gain_1)]);
% 		disp(['tot_fixed_gain_2 ' num2str(tot_fixed_gain_2)]);
% 		disp(['transmit_power ' num2str(transmit_power)]);
% 		disp(['doppler_range_correction ' num2str(doppler_range_correction)]);
% 		disp(['instrument_range_correction_tx_rx ' num2str(instrument_range_correction_tx_rx)]);
% 		disp(['instrument_range_correction ' num2str(instrument_range_correction)]);
% 		disp(['instrument_sigma0_correction_tx_rx ' num2str(instrument_sigma0_correction_tx_rx)]);
% 		disp(['instrument_sigma0_correction_rx ' num2str(instrument_sigma0_correction_rx)]);
% 		disp(['internal_phase_correction ' num2str(internal_phase_correction)]);
% 		disp(['external_phase_correction ' num2str(external_phase_correction)]);
% 		disp(['noise_power ' num2str(noise_power)]);
        
    end
    
    if(end_of_file==0)
        %-----------------------------%
        %--  read Corrections Group --%
        %-----------------------------%
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

        %-------------------------%
        %-- read waveform group --%
        %-------------------------%

        %-------------------------------------------------------%
        %-- !!! since data is no needed, we don't read it !!!!--%
        %-------------------------------------------------------%
    if(strcmp(mode,'SIN'))
        fread(fid,2621520,'int8');
    elseif(strcmp(mode,'SAR'))
        fread(fid,327760,'int8');
    end
    

        %--------------%
        %--  OUTPUT  --%
        %--------------%

        fbr_dsr = struct(   'time_orbit_group', time_orbit_group,...
                            'measurements_group',measurements_group,...
                            'corrections_group', corrections_group...
                        );

        clear time_orbit_group
        clear measurements_group
        clear corrections_group

        fbr_ds(record_num) = fbr_dsr;
    end
end
fclose(fid);


