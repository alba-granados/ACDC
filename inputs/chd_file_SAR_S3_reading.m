CHD_filename='S3A_SR_CCDB_CHAR_NOM.20151125162650_1.nc';

ncread(CHD_filename,'/SRAL_Datation_features/Time_Shift_LRM')
ncread(CHD_filename,'/SRAL_Datation_features/Time_Shift_SAR')
ncread(CHD_filename,'/SRAL_Datation_features/Nb_PRI_C_Ku')
ncread(CHD_filename,'/USO_Nominal_feature/USO_Freq_Nom')
%% SRAL_pulse_and_burst_features
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Pulse_Duration')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Sign_Slope')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/N_Burst')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Burst_Duration')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/NbPulsesBurst')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Freq_Ku')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Freq_C')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Band_Ku')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Band_C')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Abs_Ref_Track_Ku')
ncread(CHD_filename,'/SRAL_pulse_and_burst_features/Abs_Ref_Track_C')
%% SRAL_PTR_FFT_Features
ncread(CHD_filename,'/SRAL_PTR_FFT_Features/Ratio_PTR_FFT')
ncread(CHD_filename,'/SRAL_PTR_FFT_Features/FFT_Step_Freq_Ku')
ncread(CHD_filename,'/SRAL_PTR_FFT_Features/FFT_Step_Freq_C')
ncread(CHD_filename,'/SRAL_PTR_FFT_Features/FFT_Step_Time_Ku')
ncread(CHD_filename,'/SRAL_PTR_FFT_Features/FFT_Step_Time_C')
ncread(CHD_filename,'/SRAL_PTR_FFT_Features/Ind_0_Freq')
%% SRAL_Antenna
ncread(CHD_filename,'/SRAL_Antenna/Ant_3dB_BW_Ku')
ncread(CHD_filename,'/SRAL_Antenna/Ant_3dB_BW_C')
ncread(CHD_filename,'/SRAL_Antenna/Ant_Gain_Ku')
ncread(CHD_filename,'/SRAL_Antenna/Ant_Gain_C')
%% SRAL_AGC_Correction_Table
ncread(CHD_filename,'/SRAL_AGC_Correction_Table/AGC_Table_Ku')
ncread(CHD_filename,'/SRAL_AGC_Correction_Table/AGC_Table_C')
%% SRAL_Calibration_Path
ncread(CHD_filename,'/SRAL_Calibration_Path/Ratio_TRC_Ku')
ncread(CHD_filename,'/SRAL_Calibration_Path/Ratio_TRC_C')
ncread(CHD_filename,'/SRAL_Calibration_Path/Dist_Dup_Ant_Ku')
ncread(CHD_filename,'/SRAL_Calibration_Path/Dist_Dup_Ant_C')
%% SRAL_prelaunch_PTR_LRM
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/PTR_Ref_Power_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/PTR_Ref_Power_LRM_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/Half_Width_PTR_Ref_Power_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/Half_Width_PTR_Ref_Power_LRM_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/Ref_Width_Main_Lobe_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/Ref_Width_Main_Lobe_LRM_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/Ref_Diff_Trav_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/Ref_Diff_Trav_LRM_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/AGC1_PTR_PL_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/AGC1_PTR_PL_LRM_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/AGC2_PTR_PL_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/Np_PTR_PL_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/Np_PTR_PL_LRM_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/FFT_PTR_F_PL_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/FFT_PTR_F_PL_LRM_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/FFT_PTR_F_PL_LRM_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_LRM/FFT_PTR_F_PL_LRM_C')
%% SRAL_prelaunch_PTR_SAR
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/PTR_Ref_Power_SAR_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/PTR_Ref_Power_SAR_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/Half_Width_PTR_Ref_Power_SAR_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/Half_Width_PTR_Ref_Power_SAR_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/Ref_Width_Main_Lobe_SAR_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/Ref_Width_Main_Lobe_SAR_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/Ref_Diff_Trav_SAR_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/Ref_Diff_Trav_SAR_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/AGC1_PTR_PL_SAR_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/AGC1_PTR_PL_SAR_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/AGC2_PTR_PL_SAR_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/AGC2_PTR_PL_SAR_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/Np_PTR_PL_SAR_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/Np_PTR_PL_SAR_C')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/FFT_PTR_F_PL_SAR_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_PTR_SAR/FFT_PTR_F_PL_SAR_C')
%% SRAL_prelaunch_GPRW
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Width_Half_Band_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Width_Half_Band_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Mean_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Mean_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Mean_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Mean_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Std_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Std_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Std_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Std_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Diff_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Diff_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Diff_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Diff_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Slope_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Slope_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Slope_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Slope_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Std_Slope_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Std_Slope_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Std_Slope_Ku')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Right_Std_Slope_C')
ncread(CHD_filename,'/SRAL_prelaunch_GPRW/Pre_Left_Slope_C')
%% SRAL_CAL1_LRM_routine_identifiers
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Cal1_LRM_Nom_Att_Status')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Cal1_LRM_Nom_Instr_Mode')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Cal1_LRM_Nom_Nb_Cycle')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Cal1_LRM_Nom_AGC1_Ku')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Cal1_LRM_Nom_AGC1_C')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Cal1_LRM_Nom_AGC2_Ku')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Cal1_LRM_Nom_AGC2_C')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Cal1_LRM_Nom_Synth_Freq')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Np_PTR_LRM_Ku')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/Np_PTR_LRM_C')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/FFT_PTR_F_LRM_Ku')
ncread(CHD_filename,'/SRAL_CAL1_LRM_routine_identifiers/FFT_PTR_F_LRM_C')
%% SRAL_CAL1_SAR_routine_identifiers
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Cal1_SAR_Nom_Att_Status')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Cal1_SAR_Nom_Instr_Mode')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Cal1_SAR_Nom_Nb_Cycle')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Cal1_SAR_Nom_AGC1_Ku')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Cal1_SAR_Nom_AGC1_C')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Cal1_SAR_Nom_AGC2_Ku')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Cal1_SAR_Nom_AGC2_C')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Cal1_SAR_Nom_Synth_Freq')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Np_PTR_SAR_Ku')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/Np_PTR_SAR_C')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/FFT_PTR_F_SAR_Ku')
ncread(CHD_filename,'/SRAL_CAL1_SAR_routine_identifiers/FFT_PTR_F_SAR_C')
%% SRAL_CAL2_routine_identifiers
ncread(CHD_filename,'/SRAL_CAL2_routine_identifiers/Cal2_Nom_Nb_Cycle_Ku')
ncread(CHD_filename,'/SRAL_CAL2_routine_identifiers/Cal2_Nom_Nb_Cycle_C')
ncread(CHD_filename,'/SRAL_CAL2_routine_identifiers/Cal2_Nom_AGC1_Ku')
ncread(CHD_filename,'/SRAL_CAL2_routine_identifiers/Cal2_Nom_AGC1_C')
ncread(CHD_filename,'/SRAL_CAL2_routine_identifiers/Cal2_Nom_AGC2_Ku')
ncread(CHD_filename,'/SRAL_CAL2_routine_identifiers/Cal2_Nom_AGC2_C')
%% SRAL_reference_GPRW
ncread(CHD_filename,'/SRAL_reference_GPRW/GPRW_Over_Ref_Ku')
ncread(CHD_filename,'/SRAL_reference_GPRW/GPRW_Over_Ref_C')
%% SRAL_CAL1_autocal_parameters
ncread(CHD_filename,'/SRAL_CAL1_autocal_parameters/NB_AGC_Couples')
ncread(CHD_filename,'/SRAL_CAL1_autocal_parameters/Cal1_SAR_Nom_Instr_Mode')
ncread(CHD_filename,'/SRAL_CAL1_autocal_parameters/Cal1_SAR_Nom_Nb_Cycle')
ncread(CHD_filename,'/SRAL_CAL1_autocal_parameters/Cal1_SAR_Nom_Synth_Freq')
ncread(CHD_filename,'/SRAL_CAL1_autocal_parameters/Np_PTR_AGC_Ku')
ncread(CHD_filename,'/SRAL_CAL1_autocal_parameters/FFT_PTR_F_AGC_Ku')
%% SRAL_radar_data_base
dim_AGC_table_Ku  = 76;
dim_AGC_ref_table = 63;
ncread(CHD_filename,'/SRAL_radar_data_base/PRF_LRM')
ncread(CHD_filename,'/SRAL_radar_data_base/PRF_SAR')
ncread(CHD_filename,'/SRAL_radar_data_base/BRF')
ncread(CHD_filename,'/SRAL_radar_data_base/N_Pulse')
ncread(CHD_filename,'/SRAL_radar_data_base/N_Pulse_Ku')
ncread(CHD_filename,'/SRAL_radar_data_base/N_Pulse_C')
ncread(CHD_filename,'/SRAL_radar_data_base/G_CFA_Trk_LRM_Ku')
ncread(CHD_filename,'/SRAL_radar_data_base/G_CFA_Trk_LRM_C')
ncread(CHD_filename,'/SRAL_radar_data_base/G_ACCU_Trk_LRM_Ku')
ncread(CHD_filename,'/SRAL_radar_data_base/G_ACCU_Trk_LRM_C')
ncread(CHD_filename,'/SRAL_radar_data_base/G_CFA_Trk_SAR_Ku')
ncread(CHD_filename,'/SRAL_radar_data_base/G_CFA_Trk_SAR_C')
ncread(CHD_filename,'/SRAL_radar_data_base/G_ACCU_Trk_SAR_Ku')
ncread(CHD_filename,'/SRAL_radar_data_base/G_ACCU_Trk_SAR_C')
ncread(CHD_filename,'/SRAL_radar_data_base/G_CFA_Cal1_LRM_Ku')
ncread(CHD_filename,'/SRAL_radar_data_base/G_CFA_Cal1_LRM_C')
ncread(CHD_filename,'/SRAL_radar_data_base/G_ACCU_Cal1_LRM_Ku')
ncread(CHD_filename,'/SRAL_radar_data_base/G_ACCU_Cal1_LRM_C')
ncread(CHD_filename,'/SRAL_radar_data_base/Cal1_AGC_sequence')
ncread(CHD_filename,'/SRAL_radar_data_base/AGC_ref_table')
ncread(CHD_filename,'/SRAL_radar_data_base/Amb_Order')















