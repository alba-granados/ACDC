%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the UNPACKING & DECODING
% algorithm as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: Read ISP data
% 
% INPUTs : Filename to read
% OUTPUTs: TM Structure as defined on isardSAT_JasonCS_DPM
%
% ----------------------------------------------------------
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/09/2013)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Process_ID,   Process_ID2,        burst_sar_isp,              source_seq_count_sar_isp,   ...
         time_sar_isp,  inst_id_sar_isp,    trk_config_sar_isp,         tm_mode_id_sar_isp,...
         nimp_sar_isp,  pri_sar_isp,        ambiguity_order_sar_isp,...
         h0_sar_isp,    cor2_sar_isp,       loss_track_criterion_isp,...
         att_sar_ku_isp,att_sar_c_isp,      nav_bulletin_status_isp,...
         delta_alt_isp, rec_counter_sar_isp,wfm_iq_sar_ku_isp]...
            = read_ISP(...
                filename_ISP, size_ISP)       

global N_total_bursts_sar_ku_isp pri_T0_unit_conv_chd T0_nom
global N_bursts_cycle_chd sec_in_day_cst pulse_length_chd
global N_ku_pulses_burst_chd N_pulses_burst N_pulses_rc N_samples_sar_chd
global burst_duration_sar_chd inputPath N_samples_rmc_chd
global bri_nom brf_nom prf_sar_nom pri_sar_nom i_ISP


%% TM Type identification
fid_ISP             = fopen([inputPath filename_ISP],'r','b');           %opening in 'r' read mode,  'b' big endian
ID                  = fread(fid_ISP,2,'ubit8');               %[0:2 Version Number, 3 Type, 4 Data Field Header Flag, 5:7 MSB Process Identifier][0:3 LSB Process Identifier, 4:7 Packet Category]
                    frewind(fid_ISP);
ID1                 = dec2bin(ID(1),8); 
ID2                 = dec2bin(ID(2),8);
Process_ID          = bin2dec([ID1(6:8) ID2(1:4)]);
clear ID ID1 ID2

switch(Process_ID)
    case 58
        type = 1; N_total_bursts_sar_ku_isp(i_ISP)         = size_ISP/33346;       % 33346 Bytes every TM_ECHO_SAR
        
    case 59
        type = 2; N_total_bursts_sar_ku_isp(i_ISP)          = size_ISP/16454;       % 16454 Bytes every TM_ECHO_RMC
    otherwise
        disp('TM type erroneus')
end

wfm_iq_sar_ku_isp_i                 = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_ku_pulses_burst_chd,N_samples_sar_chd);
wfm_iq_sar_ku_isp_q                 = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_ku_pulses_burst_chd,N_samples_sar_chd);

% progressbar(['Reading ISP ' num2str(i_ISP) '/' num2str(size(filename_ISP,1))],[],[],[],[],[],[],[],[],[],[],[],[],[]);
for i_burst = 1:N_total_bursts_sar_ku_isp(i_ISP)
    progressbar((i_burst-0.01)/(N_total_bursts_sar_ku_isp(i_ISP)),[],[],[],[],[],[],[],[],[],[],[],[],[]);
    check = fread(fid_ISP,1,'int32');
    fseek(fid_ISP,-4,'cof');
    %------------check end of file---------------
    if(isempty(check) == 1) %If we don't find anything inside means that the file has ended
        disp('Unexpected EOF')
    else
    %         fread(fid_ISP,2,'char');
        %---UNPACKING Packet Headers Field (Field 1)
        ID                                 = fread(fid_ISP,2,'ubit8');  %[0:2 Version Number, 3 Type, 4 Data Field Header Flag, 5:7 MSB Process Identifier][0:3 LSB Process Identifier, 4:7 Packet Category]
        ID1                                = dec2bin(ID(1),8); 
        ID2                                = dec2bin(ID(2),8);
        rec_counter_sar_isp(i_burst)       = i_burst;

        %------------check TM Type---------------
%         if (bin2dec([ID1(1:3) ID2(5:7)]~=Process_ID))
%         if(iChain<3)
            Process_ID2(i_burst)                     = bin2dec([ID1(6:8) ID2(1:4)]);
%         else
%             Process_ID2(i_burst) = bin2dec(ID2);
%         end
        if (Process_ID2(1)~=Process_ID)
            disp('TM type different from first TM') 
        else
            clear ID ID1 ID2
            source_seq_count_sar_isp(i_burst)      = fread(fid_ISP,1,'uint16');             %[0:1 Grouping Flag, 2:7 MSB sequence count][0:7 LSB Sequence Count]
            packet_length(i_burst)                 = fread(fid_ISP,1,'uint16');             %Number of bytes of the TM (except packet header) minus 1 Byte 

            %% UNPACKING Data Field Header (Field 2)

            data_field_header                      = fread(fid_ISP,4,'char');             %[All the header]

            day                                         = fread(fid_ISP,1,'int16');            % Day Epoch starting at 1.1.2000 
            milliseconds                                = fread(fid_ISP,1,'int32');            % Milliseconds of the day
            microseconds                                = fread(fid_ISP,1,'int16');            % Microseconds of milliseconds
            %% DECODING Time (Field 2)
            time_sar_isp(i_burst) = day * sec_in_day_cst + milliseconds * 1e-3 + microseconds * 1e-6 ;
            %% UNPACKING Application Data Common Part  (Fields 3-22)
                                                          fseek(fid_ISP,1,'cof');              % Field 3 Reserved 
            inst_id_sar_isp(i_burst)                    = fread(fid_ISP,1,'unsigned char');    % Field 4 P4 identifier [0:6 Reserved, 7 P4 ID] 0= nominal, 1=redundant
            trk_config_sar_isp(i_burst)                 = fread(fid_ISP,1,'unsigned char');    % Field 5 Tracking Configuration
            tm_mode_id_sar_isp(i_burst)            = fread(fid_ISP,1,'unsigned char');    % Field 6 TM Mode identifier

            nimp_sar_isp(i_burst)                  = fread(fid_ISP,1,'ushort');           % Field 7 Number of pulses
            pri_sar_isp(i_burst)                   = fread(fid_ISP,1,'ushort');           % Field 8 Pulse Repetition Interval (T0*8 units)
            ambiguity_order_sar_isp(i_burst)       = fread(fid_ISP,1,'ushort');           % Field 9 Ambiguity Order (PRI units)
            h0_comp_sar_isp                             = fread(fid_ISP,1,'ulong');            % Field 10 H0 computed from navigation and DEM (T0/64 units)
            cor2_comp_sar_isp                           = fread(fid_ISP,1,'short');            % Field 11 COR2 computed from navigation and DEM (T0/64/16 units)
            h0_sar_isp(i_burst)                    = fread(fid_ISP,1,'ulong');             % Field 12 H0 applied (T0/64 units)
            cor2_sar_isp(i_burst)                  = fread(fid_ISP,1,'short');            % Field 13 COR2 applied (T0/64/16 units)
                                                          fseek(fid_ISP,1,'cof');            % Field 14 Reserved
                                                          
            loss_track_criterion_isp(i_burst)           = fread(fid_ISP,1,'char');             % Field 15 Loss of track criterion computed on the echo of the cycle (N-2) in OL mode (current cycle = N) (Bit 7): 0 = normal 1 = loss of track
            dist_err                                    = fread(fid_ISP,1,'long');            % Field 16 Distance error (T0/64 units) computed on the echo of the cycle (N-2) in OL mode (current cycle =N) 
            att_sar_ku_isp(i_burst)                = fread(fid_ISP,1,'unsigned char');    % Field 17 ATTCODE Ku band (dB units)
            att_sar_c_isp(i_burst)                 = fread(fid_ISP,1,'unsigned char');    % Field 18 ATTCODE C  band (dB units)
            ASS32                                       = fread(fid_ISP,4,'char');             % Field 19 AS32 of the previous cycle
                                                          fseek(fid_ISP,1,'cof');            % Field 20 Reserved 
            nav_bulletin_status_isp(i_burst)       = fread(fid_ISP,1,'unsigned char');    % Field 21 Navigation bulletin status (Bit 7): 0 = bulletin OK; 1 = bulletin KO
            time_nav_bulletin_isp                       = fread(fid_ISP,8,'unsigned char');    % Field 22 Time of current navigation bulletin
            
            %% UNPACKING Application Data SAR RAW      (Fields 23-26)   
            if(type==1)
                wfm_iq_sar_ku_isp_aux            = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_ku_pulses_burst_chd*N_samples_sar_chd*2);
                delta_alt_isp(i_burst)           = 0;
                                                          fseek(fid_ISP,1,'cof');             %Field 23 Reserved for SAR                            % 
                burst_sar_isp(i_burst)              	= fread(fid_ISP,1,'uint8');           % Field 24 burst number (from 1 to 7)
                wfm_iq_sar_c_isp_aux(i_burst,:)            	= fread(fid_ISP,N_samples_sar_chd*2,'int8');         % Field 25 65x256 (I,Q) samples: 1 C pulse + 64 Ku pulses [Y+V,M]
                wfm_iq_sar_ku_isp_aux(i_burst,:)            = fread(fid_ISP,N_samples_sar_chd*N_ku_pulses_burst_chd*2,'int8');      % Field 26 65x256 (I,Q) samples: 1 C pulse + 64 Ku pulses [Y+V,M]
                
                CRC                                         = fread(fid_ISP,1,'uint16');             % Field 26 CRC
                for i_pulse = 1:N_ku_pulses_burst_chd
                    for i_sample = 1:N_samples_sar_chd
%                         wfm_iq_sar_ku_isp_i(i_burst,i_pulse,i_sample) = wfm_iq_sar_ku_isp_aux(i_burst,(i_pulse-1)*N_samples_sar_chd*2 + 2*i_sample-1);
%                         wfm_iq_sar_ku_isp_q(i_burst,i_pulse,i_sample) = wfm_iq_sar_ku_isp_aux(i_burst,(i_pulse-1)*N_samples_sar_chd*2 + 2*i_sample);
                        wfm_iq_sar_ku_isp_q(i_burst,i_pulse,i_sample) = wfm_iq_sar_ku_isp_aux(i_burst,(i_pulse-1)*N_samples_sar_chd*2 + 2*i_sample-1);
                        wfm_iq_sar_ku_isp_i(i_burst,i_pulse,i_sample) = wfm_iq_sar_ku_isp_aux(i_burst,(i_pulse-1)*N_samples_sar_chd*2 + 2*i_sample);
                    end
                end

            %% UNPACKING Application Data SAR RMC      (Fields 23-27)     
            elseif(type==2)
                wfm_iq_sar_ku_isp_aux            = zeros(N_total_bursts_sar_ku_isp(i_ISP),N_ku_pulses_burst_chd*N_samples_rmc_chd*2);
                
                delta_alt_isp(i_burst)                 = 1e-7*fread(fid_ISP,1,'long');             % Field 23 10-7m (from 1 to 7)
                scaling_factor_iq             = fread(fid_ISP,1,'uchar');            % Field 24 I and Q scaling factor
                burst_sar_isp(i_burst)                 = fread(fid_ISP,1,'uint8');            % Field 25 burst number (from 1 to 7)
                wfm_iq_sar_ku_isp_aux(i_burst,:)            = fread(fid_ISP,N_samples_rmc_chd*N_ku_pulses_burst_chd*2,'int8');        % Field 26 64x128 (I,Q) samples:64 Ku pulses [Y,R]
                CRC                                         = fread(fid_ISP,1,'uint16');             % Field 27 CRC
                for i_pulse = 1:N_ku_pulses_burst_chd
                    for i_sample = 1:N_samples_rmc_chd
%                         wfm_iq_sar_ku_isp_i(i_burst,i_pulse,i_sample) = scaling_factor_iq(i_burst)*wfm_iq_sar_ku_isp_aux(i_burst,(i_pulse-1)*N_samples_sar_chd*2 + 2*i_sample-1);
%                         wfm_iq_sar_ku_isp_q(i_burst,i_pulse,i_sample) = scaling_factor_iq(i_burst)*wfm_iq_sar_ku_isp_aux(i_burst,(i_pulse-1)*N_samples_sar_chd*2 + 2*i_sample);
                        wfm_iq_sar_ku_isp_q(i_burst,i_pulse,i_sample) = scaling_factor_iq*wfm_iq_sar_ku_isp_aux(i_burst,(i_pulse-1)*N_samples_rmc_chd*2 + 2*i_sample-1);
                        wfm_iq_sar_ku_isp_i(i_burst,i_pulse,i_sample) = scaling_factor_iq*wfm_iq_sar_ku_isp_aux(i_burst,(i_pulse-1)*N_samples_rmc_chd*2 + 2*i_sample);
                    end
                end
            end 
            clear day milliseconds microseconds 
            clear h0_comp_sar_isp cor2_comp_sar_isp dist_err ASS32 time_nav_bulletin_isp CRC
            
        end
    end
end

nimp_sar_isp = zeros(1,N_total_bursts_sar_ku_isp(i_ISP)) + 455; %RMC
%Once the RMC samples have been read from the ISP, N = 128 is not useful anymore
% N_samples_sar_chd = N_samples_sar_chd;
% N_total_bursts_sar_ku_isp = N_total_bursts_sar_ku_isp;
% T0_nom = T0_nom;

N_pulses_burst = 65;


wfm_iq_sar_ku_isp = wfm_iq_sar_ku_isp_i + 1i*wfm_iq_sar_ku_isp_q;



N_pulses_rc = N_pulses_burst * N_bursts_cycle_chd;
pri_sar_nom = (pri_sar_isp * pri_T0_unit_conv_chd * T0_nom);
burst_duration_sar_chd = N_ku_pulses_burst_chd * pri_sar_nom + pulse_length_chd;
prf_sar_nom = 1./pri_sar_nom;
bri_nom = N_pulses_burst .* pri_sar_nom;
brf_nom = 1 ./ bri_nom;


fclose('all');



end



