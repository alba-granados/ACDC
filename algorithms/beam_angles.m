%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
%
% ---------------------------------------------------------
% Objective: The purpose of the beam angles is to compute, for every burst, 
% the angles between the nadir direction and the direction defined by 
% the satellite location and each surface location under the satellite's boresight. 
% 
% ----------------------------------------------------------
% Author:    Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
%
% Version  record
% 1.0 2015/01/01 Imported code from S6
% 1.1 2016/03/03 Changed pri_sar_nom for L1A.pri_sar
% 2.0 2016/04/01 Beam angles per burst
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%____________________________For each burst________________________________
function [L1A]       = beam_angles (L1A,L1BS,N_total_surf_loc,i_surf)
                                                                                            
                                                                                            
    global c_cst pi_cst 
    global freq_ku_chd N_ku_pulses_burst_chd pri_sar_nom
    global zp_fact_azimut_cnf

    first_beam_ang_index        = 1;
    
    beam_ang                    = 1./zeros(1,N_ku_pulses_burst_chd*zp_fact_azimut_cnf);
    surf_loc_index              = zeros(1,N_ku_pulses_burst_chd*zp_fact_azimut_cnf);
    
        
    vel_sat = [L1A.x_vel_sat_sar,L1A.y_vel_sat_sar,L1A.z_vel_sat_sar];

    %% 1. Computation of the boresight angles qmin and qmax
    qmin = acos (c_cst / L1A.pri_sar / (4 * freq_ku_chd * norm(vel_sat)));
    qmax = pi_cst - qmin;


    %% 2. Surface locations loop
%     i_surf = first_beam_ang_index;
    b = qmin + 1; %initialization of the variable in order to perform a do-while

    while b >= qmin && i_surf <= N_total_surf_loc

        %% 3. Beam angle computation
        v = [L1BS(i_surf).x_surf, L1BS(i_surf).y_surf, L1BS(i_surf).z_surf] - [L1A.x_sar_sat, L1A.y_sar_sat, L1A.z_sar_sat];
        w = vel_sat;

        b = acos ((v*w') / (norm(v) * norm(w)));

        %% 4. Check angular limits
        if b <= qmax && b >= qmin

            %% 5. Store the beam angle and its index
            %******************** Sanity Check **********************
            if  L1A.N_beams_sar_ku < N_ku_pulses_burst_chd*zp_fact_azimut_cnf
                L1A.N_beams_sar_ku = L1A.N_beams_sar_ku + 1;
                beam_ang(L1A.N_beams_sar_ku) = b;
                surf_loc_index(L1A.N_beams_sar_ku) = i_surf;
            else
                clear beam_ang_aux surf_loc_index_aux
                beam_ang_65 = b;                 % afegim el beam 65
                beam_ang_aux = beam_ang(2:L1A.N_beams_sar_ku);    % i ens carreguem el 1r de tots (el de m�s enrere) i aix� en tenim 64 
                beam_ang(1:L1A.N_beams_sar_ku-1) = beam_ang_aux;
                beam_ang(L1A.N_beams_sar_ku)     = beam_ang_65;
                surf_loc_index_65 = i_surf;         % fem el mateix amb els �ndexs de superf�cies
                surf_loc_index_aux = surf_loc_index(2:L1A.N_beams_sar_ku);
                surf_loc_index(1:L1A.N_beams_sar_ku-1) = surf_loc_index_aux;
                surf_loc_index(L1A.N_beams_sar_ku)     = surf_loc_index_65;
            end
        elseif b > qmax
            first_beam_ang_index = first_beam_ang_index + 1;
        end

        i_surf = i_surf + 1;

    end
    L1A.beam_ang_nadir_index = 0;

    %order within stack to be used in the future
    % force that the central beam to the closest direction perpendicular to the velocity and forward (negative).
    L1A.beam_ang_index= N_ku_pulses_burst_chd/2*zp_fact_azimut_cnf+1;
%     [~,beams_shift] = min(abs((beam_ang(:)-pi_cst/2)));
%     if((beam_ang(beams_shift)-pi_cst/2)>0)
%         beams_shift = beams_shift+1;
%     end
%     beams_shift = N_ku_pulses_burst_chd/2+1-beams_shift;
%     


%     L1A.start_beam     = max(1,beams_shift+1);
%     L1A.end_beam       = min(beams_shift+ L1A.N_beams_sar_ku,N_ku_pulses_burst_chd);
%     L1A.beam_index  = L1A.start_beam:L1A.end_beam;
%         L1A.beam_ang(L1A.beam_index) = beam_ang(1:L1A.N_beams_sar_ku);
%         L1A.surf_loc_index(L1A.beam_index) = surf_loc_index(1:L1A.N_beams_sar_ku);
    L1A.beam_ang=beam_ang;
    L1A.surf_loc_index=surf_loc_index;

    
end



