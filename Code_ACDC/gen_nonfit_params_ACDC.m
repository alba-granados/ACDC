function nf_p = gen_nonfit_params_ACDC (L1BS,L1B) 
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This function initializates the non fitting parameters' structure
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT v1 29/09/2016
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       data    =   input data structure for the L2 processor
% OUTPUT:
%       nf_p           =   Non-fitting params

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - 
% Versions control:
% v1.0:
%% -------------- Global variables ----------------------------------------
% ---------------- Constants ----------------------------------------------
global c_cst
global semi_major_axis_cst semi_minor_axis_cst
% ---------------- Configuration parametres -------------------------------
global antenna_beamwidth_alt_ku_chd antenna_beamwidth_act_ku_chd
global freq_ku_chd hamming_window_cnf window_rg_cnf
global N_ku_pulses_sar_chd N_bursts_cycle_chd bw_ku_chd N_samples

%-------------------- Configuration of retracker ACDC ---------------------
global mission cnf_p_ACDC

%% ------------- Computation of the burst nadir/related parameters --------
%Compute nadir burst for a given surface
[~,beam_index_nadir]=min(abs(-pi/2+L1BS.beam_ang_surf'));
%burst_index_nadir= L1BS.burst_index(beam_index_nadir);
nf_p.pri_surf            = L1BS.pri_sar_sat_beam(beam_index_nadir);


%% ------------------- Height statistical distribution --------------------
%if ~cnf_p.rou_flag
    nf_p.waveskew   =   0;
    nf_p.EMbias     =   0;
    nf_p.rou    =   cnf_p_ACDC.rou; % from Yuguan Liu
%end

%% ------------------ Array values ----------------------------------------
%% ------------------ Geometry/Projections --------------------------------
% Height
nf_p.h       =   L1BS.alt_sat;
% Velocity
nf_p.v_sat   =   norm([L1BS.x_vel_sat,L1BS.y_vel_sat,L1BS.z_vel_sat]);
% Spherical earth term
R_Earth =sqrt((semi_major_axis_cst^4*cos(L1BS.lat_sat*pi/180).^2+semi_minor_axis_cst^4*sin(L1BS.lat_sat*pi/180).^2)./((semi_major_axis_cst*cos(L1BS.lat_sat*pi/180)).^2+(semi_minor_axis_cst*sin(L1BS.lat_sat*pi/180)).^2));
alpha        =   1 + (nf_p.h)./(R_Earth);

% Pitch and roll angles projection as distances
% assuming no change of sign is performed in the lecture
%take into account that the reference system for Model (x,y,z) is having z
%normal to the surface and x in the plane formed by satellite velocity and
%z : positive y component on the left_side; x positive forward movement
%satellite
switch mission
    case {'CS2','CR2'}
        % Satellite-frame coordinates: x-along-track, z-opposite direction
        % to nadir and positive y on left-side of satellite
        %pitch is defined as positive/negative nose down/up respectively        
        nf_p.xp      =   L1BS.alt_sat.*(-1.0*L1BS.pitch_surf);
        %roll-angle is defined as positive/negative antenna right down/up        
        nf_p.yp      =   L1BS.alt_sat.*(1.0*L1BS.roll_surf);
    case {'S3','S3_'}  % alba: I followed SAR altimeter analytical open ocean conventional retracker/algorithms/gen_nonfit_params_EM.m
        %TBD
        %pitch is defined as positive/negative nose up/down respectively (clockwise w.r.t Y-axis looking from origin )
        nf_p.xp      =   L1BS.alt_sat.*(1.0*L1BS.pitch_surf);
        %roll-angle is defined as positive/negative satellite-right-side
        %down/up (clockwise angle w.r.t x-axis looking from origin)
        nf_p.yp      =   L1BS.alt_sat.*(1.0*L1BS.roll_surf);
    case {'S6','S6_'}
        % Satellite-frame coordinates: x-along-track, z direction
        % to nadir and positive y on right-side of satellite
        %pitch is defined as positive/negative nose up/down respectively (clockwise w.r.t Y-axis looking from origin )
        nf_p.xp      =   L1BS.alt_sat.*(1.0*L1BS.pitch_surf);
        %roll-angle is defined as positive/negative satellite-right-side
        %down/up (clockwise angle w.r.t x-axis looking from origin)
        nf_p.yp      =   L1BS.alt_sat.*(L1BS.roll_surf);        
end


%Equivalent antenna beamwidth projection on ground
nf_p.alphax  =   (8*log(2))./(L1BS.alt_sat.*antenna_beamwidth_alt_ku_chd).^2;
nf_p.alphay  =   (8*log(2))./(L1BS.alt_sat.*antenna_beamwidth_act_ku_chd).^2; %8 comes from fact that using the two way antenna pattern

%% --------------- Sampling spacing ---------------------------------------
% along-track
nf_p.Lx      =   c_cst*L1BS.alt_sat./(2*nf_p.v_sat.*freq_ku_chd*N_ku_pulses_sar_chd.*nf_p.pri_surf);
% across-track
%nf_p.Ly      =   sqrt(c_cst*data.GEO.H./(alpha*bw_rx_ku_chd*cnf_p.ZP));
nf_p.Ly      =   sqrt(c_cst*L1BS.alt_sat./(alpha*bw_ku_chd));
% vertical
%nf_p.Lz      =   0.5*c_cst/(bw_rx_ku_chd*cnf_p.ZP);
nf_p.Lz      =   0.5*c_cst/(bw_ku_chd);

%nf_p.oversampling=fs_clock_ku_chd/bw_ku_chd; %due to differences in BW and fs (like in Sentinel-6)

% not clear from where this term is coming
nf_p.Lgamma  =   alpha./(2*L1BS.alt_sat.*nf_p.alphay);

% %% --------------- Number of looks that form the stack --------------------
% % it would be probably more convenient to define a unique structure value of Neff
% % no need to use data.SIN.Neff (will be kept as originally proposed just in case there are some future elaborations on the SARin data)
nf_p.Neff=L1B.N_beams_start_stop;


%% ------------------- PTR gaussian approximation -------------------------
%first approach assuming no hamming in range: boxcar (to be updated)
if window_rg_cnf
    A_s2Gr_chd=1.0101;
    alpha_gr_chd=1.0/(2*(0.59824).^2);  
    %alpha_gr_chd=1.0/(2*(0.5408).^2);  
else
    A_s2Gr_chd=1.0196;
    alpha_gr_chd=1.0/(2*(0.36012).^2);    
end

% if applied right now only on azimuth
if hamming_window_cnf
    A_s2Ga_chd=1.0101;
    alpha_ga_chd=1.0/(2*(0.59824).^2);    
    %alpha_ga_chd=1.0/(2*(0.5408).^2);    
else
    A_s2Ga_chd=1.0196;
    alpha_ga_chd=1.0/(2*(0.36012).^2);
end

%% ------------------- Constant Values ------------------------------------
nf_p.Npulses    =   N_ku_pulses_sar_chd;
nf_p.alphag_a   =   alpha_ga_chd;
nf_p.alphag_r   =   alpha_gr_chd;
nf_p.A_s2Ga     =   A_s2Ga_chd;
nf_p.A_s2Gr     =   A_s2Gr_chd;
nf_p.Nbcycle    =   N_bursts_cycle_chd;
nf_p.bw_Rx      =   bw_ku_chd;
nf_p.Nsamples   =   N_samples;



end