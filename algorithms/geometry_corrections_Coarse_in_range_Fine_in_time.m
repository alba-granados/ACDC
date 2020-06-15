%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the GEOMETRY CORRECTIONS as described in the
% isardSAT_JasonCS_DPM_JC-DS-ISR-SY-0006_v5a_20130605
%
% ---------------------------------------------------------
% Objective: The purpose of the range compression is to perform range 
% compression of the input bursts (with an FFT) and then generate the 
% power waveforms.
% 
% INPUTs: 
%     Name                          Units   Type                                            Origin
%     N_total_surf_loc                      do                                              Surface Locations (�7.4)

%     L1BS.beams_surf            L1BS.N_beams_stack * N_samples_sar_chd * do         Azimuth Processing (�7.9)
%     L1BS.N_beams_stack                        do                                              Beam Angles (�7.8)
%     N_samples_sar_chd                     us                                              CHD file (�2.2.1)
%     L1BS.x_surf                        m       do                                              Surface Locations (�7.4)
%     L1BS.y_surf                        m       do                                              Surface Locations (�7.4)
%     L1BS.z_surf                        m       do                                              Surface Locations (�7.4)





%     L1A.x_sar_sat                     m       N_total_bursts_sar_ku * do                      Surface Locations (�7.4)
%     L1A.y_sar_sat                     m       N_total_bursts_sar_ku * do                      Surface Locations (�7.4)
%     L1A.z_sar_sat                     m       N_total_bursts_sar_ku * do                      Surface Locations (�7.4)
%     L1A.x_vel_sat_sar                 m/s     N_total_bursts_sar_ku * do                      Surface Locations (�7.4)
%     L1A.y_vel_sat_sar                 m/s     N_total_bursts_sar_ku * do                      Surface Locations (�7.4)
%     L1A.z_vel_sat_sar                 m/s     N_total_bursts_sar_ku * do                      Surface Locations (�7.4)
%     L1A.win_delay_sar_ku                 s       N_total_bursts_sar_isp*N_ku_pulses_sar_chd*do   Preliminary Window Delay (�7.3)
%     L1BS.win_delay_surf                s       do                                              Surface Locations (�7.4)

%     L1BS.beam_ang_surf                 rad     L1BS.N_beams_stack * do                             Azimuth Processing (�7.9)
%     surf_loc_index                        N_beams_ sar_ku * do                            Beam Angles (�7.8)
%     L1BS.burst_index                           L1BS.N_beams_stack * do                              Azimuth Processing (�7.9)
%     pulse_length_chd              s       do                                              CHD file (�2.2.1)
%     bw_ku_chd                     Hz      do                                              CHD file (�2.2.1)
%     wv_length_ku                  m       fl                                              CHD decoding (�3.2)
%     c_cst                         m/s     do                                              CST file (�2.2.3)




% 
% OUTPUTs:  
%     Name                          Units   Type                                            Destination
%     L1BS.doppler_corr                  samples	L1BS.N_beams_stack * do                             L1B-S (�8.3.2)
%     L1BS.range_sat_surf                m       L1BS.N_beams_stack * do                             Waveforms Scaling Factor (�7.13)
%     L1BS.slant_range_corr              samples	L1BS.N_beams_stack * do                             L1B-S (�8.3.2)  
%     shift_wd_surf                 samples	L1BS.N_beams_stack * do                             L1B-S (�8.3.2) 
%     L1BS.wfm_geo_corr                  L1BS.N_beams_stack * N_samples_sar _chd * do                Range Compression (�7.11) and L1B-S (�8.3.2) 
% ----------------------------------------------------------
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (02/06/2015)
% 
% 
% The purpose of the geometry corrections is to compute and apply:
% > Doppler correction
% > Slant range correction.
% > Window delay Misalignments
%
% Version  record
% 1.0 2015/01/01    Imported code from S6
% 1.1 2015/06/01    Added computation of N_windows.
%                   Application of shif splitted in coarse and fine to avoid wrapps   
% 1.2 2016/02/18    Added zp_fact_azimut_cnf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [L1BS]  = geometry_corrections (L1A,L1BS)
                                                                                             
global pulse_length_chd N_samples_sar_chd bw_ku_chd
global wv_length_ku zp_fact_range_cnf
global c_cst pi_cst prf_sar_nom N_total_bursts_sar_ku 
global win_del_ref_chd avoid_wrapps_cnf T0_chd 


% Pre-allocating memory
doppler_range           = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));
L1BS.doppler_corr            = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));
i_samples               = 0:(N_samples_sar_chd-1);
L1BS.range_sat_surf          = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));
L1BS.slant_range_corr_time   = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));
L1BS.slant_range_corr        = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));

L1BS.wd_corr                 = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));


win_delay_ref           = zeros(1,L1BS.N_total_surf_loc);
win_delay_ref_index     = zeros(1,L1BS.N_total_surf_loc);
win_delay_ref_index_2   = zeros(1,L1BS.N_total_surf_loc);
beam_power              = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));
beam_power2             = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));
L1BS.shift                   = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));
L1BS.shift_coarse            = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack));
L1BS.good_samples            = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack),N_samples_sar_chd*zp_fact_range_cnf);
L1BS.range_sat_surf_amb      = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack)+100);



for i_burst=1:N_total_bursts_sar_ku
    vel_sat = [L1A.x_vel_sat_sar(i_burst), L1A.y_vel_sat_sar(i_burst), L1A.z_vel_sat_sar(i_burst)];
    vel_sat_norm(i_burst) = (norm(vel_sat));
end
% Doppler_bwd = 1.35/180*pi/wv_length_ku*2*max(vel_sat_norm);
% N_amb = ceil((1-min(prf_sar_nom)/Doppler_bwd )*64/2); % N_amb beams per banda.
% N_amb_stack = 7*N_amb;


 METHOD = 0; %0 = normal alignment
             %1 = optimal alignment (wrt max beam accumulated power)
             %2 = force the alignment to the closest of a given window_delay
             %3 = force to the first one.
             %4 = force to the minimum one.
             
for i_surf = 1:L1BS.N_total_surf_loc
    
    progressbar([],[],[],[],[],[],[],[],[],[],[],(i_surf)/L1BS.N_total_surf_loc,[],[],[]);
    
   
    %% 1. BEAM POWER
    if METHOD == 0
        win_delay_ref(i_surf) = L1BS.win_delay_surf(i_surf);
    elseif METHOD ==1
        beams_surf_aux1_fft = fftshift(fft(squeeze(L1BS.beams_surf(i_surf,:,:)).'),1).';
        for i_beam = 1:L1BS.N_beams_stack(i_surf)
            beam_power(i_surf,i_beam) = sum(abs(beams_surf_aux1_fft(i_beam,:)).^2);
            beam_power2(i_surf,i_beam) = max(abs(beams_surf_aux1_fft(i_beam,:)).^2);
        end
        [~,win_delay_ref_index(i_surf)] = max(beam_power(i_surf,1:L1BS.N_beams_stack(i_surf)));
        win_delay_ref(i_surf) = L1A.win_delay_sar_ku(L1BS.burst_index(i_surf,win_delay_ref_index(i_surf)));
    elseif METHOD==2
        
        win_delay_ref(i_surf) = win_del_ref_chd;
        [~,win_delay_ref_index(i_surf)] = min(abs(win_delay_ref(i_surf)-L1A.win_delay_sar_ku(L1BS.burst_index(i_surf,1:L1BS.N_beams_stack(i_surf)))));
    elseif METHOD==3 
        win_delay_ref_index(i_surf)=1;
        win_delay_ref(i_surf) = L1A.win_delay_sar_ku(L1BS.burst_index(i_surf,win_delay_ref_index(i_surf)));
    elseif METHOD==4
        [win_delay_ref(i_surf),win_delay_ref_index(i_surf)] = min(L1A.win_delay_sar_ku(L1BS.burst_index(i_surf,1:L1BS.N_beams_stack(i_surf))));
                       
    end
      
    for i_beam = 1:L1BS.N_beams_stack(i_surf)
        %% 2. DOPPLER CORRECTION
        vel_sat  = [L1A.x_vel_sat_sar(L1BS.burst_index(i_surf,i_beam)), L1A.y_vel_sat_sar(L1BS.burst_index(i_surf,i_beam)), L1A.z_vel_sat_sar(L1BS.burst_index(i_surf,i_beam))];
        % Range correction computation
        doppler_range(i_surf,i_beam) = (- c_cst / wv_length_ku * norm(vel_sat) * cos(L1BS.beam_ang_surf(i_surf,i_beam))) * pulse_length_chd / bw_ku_chd;
        % Doppler correction computation
        L1BS.doppler_corr(i_surf,i_beam) = 2 / c_cst * doppler_range(i_surf,i_beam) / L1BS.T0_sar_surf(i_surf,i_beam);
%         L1BS.doppler_corr(i_surf,i_beam) = 0; %testing
        
        %% 3. SLANT RANGE CORRECTION
        % Range computation
        L1BS.range_sat_surf(i_surf,i_beam) = norm ([L1A.x_sar_sat(L1BS.burst_index(i_surf,i_beam)), ...
                                               L1A.y_sar_sat(L1BS.burst_index(i_surf,i_beam)), ...
                                               L1A.z_sar_sat(L1BS.burst_index(i_surf,i_beam))] - [L1BS.x_surf(i_surf), L1BS.y_surf(i_surf), L1BS.z_surf(i_surf)]);
        % Range delay computation
        L1BS.slant_range_corr_time(i_surf,i_beam) = win_delay_ref(i_surf) - (L1BS.range_sat_surf(i_surf,i_beam) * 2 / c_cst);
        L1BS.slant_range_corr(i_surf,i_beam) = L1BS.slant_range_corr_time(i_surf,i_beam) / L1BS.T0_sar_surf(i_surf,i_beam);
    
        %% 4. WINDOW DELAY MISALIGNMENTS CORRECTION
        L1BS.wd_corr(i_surf,i_beam) = - (win_delay_ref(i_surf) - L1A.win_delay_sar_ku(L1BS.burst_index(i_surf,i_beam))) / L1BS.T0_sar_surf(i_surf,i_beam);
    end
    
    L1BS.shift(i_surf,:) = L1BS.doppler_corr(i_surf,:) + ...
                           L1BS.slant_range_corr(i_surf,:) + ...
                           L1BS.wd_corr(i_surf,:);
    [~,win_delay_ref_index_2(i_surf)] = min(L1BS.shift(i_surf,:));
    L1BS.shift(i_surf,:) = L1BS.shift(i_surf,:)- L1BS.shift(i_surf,win_delay_ref_index_2(i_surf));

end
L1BS.shift_coarse = round(L1BS.shift);
L1BS.shift_fine = L1BS.shift-L1BS.shift_coarse;


%% increase the window size 
if(avoid_wrapps_cnf)
    L1BS.N_windows  = ceil(max(max(abs(L1BS.shift))/N_samples_sar_chd)+1);
    truncated_zeros = zeros(1,N_samples_sar_chd*zp_fact_range_cnf*(L1BS.N_windows-1));
else
    L1BS.N_windows  = 1;
    truncated_zeros =[];
end

L1BS.beam_geo_corr          = zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack),N_samples_sar_chd*zp_fact_range_cnf*L1BS.N_windows);

 if(strcmp(L1A.mode,'SIN'))
    L1BS.beam_geo_corr_2 	= zeros(L1BS.N_total_surf_loc,max(L1BS.N_beams_stack),N_samples_sar_chd*zp_fact_range_cnf*L1BS.N_windows);
 end



%% 5. APPLYING THE FINE CORRECTIONS IN TIME 
for i_surf = 1:L1BS.N_total_surf_loc
    progressbar([],[],[],[],[],[],[],[],[],[],[],(i_surf)/L1BS.N_total_surf_loc,[],[],[]);
    
    wfm_geo_corr_aux = zeros (max(L1BS.N_beams_stack),N_samples_sar_chd);
    if(strcmp(L1A.mode,'SIN'))
        wfm_geo_corr_aux_2 = zeros (max(L1BS.N_beams_stack),N_samples_sar_chd);
    end
    
    
    for i_beam = 1:L1BS.N_beams_stack(i_surf)
      
        wfm_geo_corr_aux(i_beam,:) = squeeze(L1BS.beams_surf(i_surf,i_beam,:)).' .* exp(2i*pi_cst/N_samples_sar_chd*L1BS.shift_fine(i_surf,i_beam).*i_samples);
        beam_zp_fft                     = fftshift(fft(wfm_geo_corr_aux(i_beam,:).', N_samples_sar_chd * zp_fact_range_cnf),1).'/sqrt(N_samples_sar_chd);
        beam_surf_aux1_fft_increased    = (cat(2,beam_zp_fft,truncated_zeros));
        beam_surf_aux1_fft_aligned      = circshift(beam_surf_aux1_fft_increased,[0,L1BS.shift_coarse(i_surf,i_beam)*zp_fact_range_cnf]);
        L1BS.beam_geo_corr(i_surf,i_beam,:) = beam_surf_aux1_fft_aligned;
        
        if(strcmp(L1A.mode,'SIN'))
         	wfm_geo_corr_aux_2(i_beam,:)    = squeeze(L1BS.beams_surf_2(i_surf,i_beam,:)).' .* exp(2i*pi_cst/N_samples_sar_chd*L1BS.shift_fine(i_surf,i_beam).*i_samples);
            beam_zp_fft                     = fftshift(fft(wfm_geo_corr_aux_2(i_beam,:).', N_samples_sar_chd * zp_fact_range_cnf),1).'/sqrt(N_samples_sar_chd);
            beam_surf_aux1_fft_increased    = (cat(2,beam_zp_fft,truncated_zeros));
            beam_surf_aux1_fft_aligned      = circshift(beam_surf_aux1_fft_increased,[0,L1BS.shift_coarse(i_surf,i_beam)*zp_fact_range_cnf]);
            L1BS.beam_geo_corr_2(i_surf,i_beam,:) = beam_surf_aux1_fft_aligned;
        end
        
        
%        if abs(L1BS.shift_coarse(i_surf,i_beam))* zp_fact_range_cnf >= N_samples_sar_chd* zp_fact_range_cnf
%          L1BS.good_samples(i_surf,i_beam,:)=0;
%        else
        if(L1BS.shift_coarse(i_surf,i_beam)>0)
            start_sample=ceil(L1BS.shift_coarse(i_surf,i_beam))*zp_fact_range_cnf+1;
            final_sample=N_samples_sar_chd*zp_fact_range_cnf;
        else
            start_sample=1;
            final_sample=(N_samples_sar_chd*zp_fact_range_cnf+floor((L1BS.shift_coarse(i_surf,i_beam))*zp_fact_range_cnf));
        end

        L1BS.good_samples(i_surf,i_beam,start_sample:final_sample)=1;
        
        
    end
    
    
   
end


% L1BS = rmfield(L1BS,'beams_surf');
if(strcmp(L1A.mode,'SIN'))
%     L1BS = rmfield(L1BS,'beams_surf_2');
end
end

    


