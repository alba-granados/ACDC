%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SURFACE & WINDOW DELAY ALGINMENT %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of the surface & window delay alignment is to coregister the 
% different leading edge positions of the surfaces over the track. Need to
% align also the window delay.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%____________________________For each stack________________________________

function [L1BS,L1B] = surface_win_delay_alignment (L1BS,L1B,win_delay_surf_ref,alt_sat_ref)

global zp_fact_range_cnf T0_chd c_cst mode processing_mode_cnf

wd_shift = ((L1BS.win_delay_surf-win_delay_surf_ref)-(L1BS.alt_sat-alt_sat_ref)*2/c_cst) / T0_chd;
L1B.wfm_cor_i2q2_sar_ku_wdcorr = circshift(L1B.wfm_cor_i2q2_sar_ku,[0,round(wd_shift*zp_fact_range_cnf)]);
L1BS.win_delay_surf_aligned = L1BS.win_delay_surf-wd_shift*T0_chd;
if strcmp(mode,'SIN') && strcmp(processing_mode_cnf,'SIN')
   L1B.wfm_cor_i2q2_sar_ku_2_wdcorr = circshift(L1B.wfm_cor_i2q2_sar_ku_2,[0,round(wd_shift*zp_fact_range_cnf)]);
   L1B.phase_difference_wdcorr = circshift(L1B.phase_difference,[0,round(wd_shift*zp_fact_range_cnf)]);
   L1B.coherence_wdcorr = circshift(L1B.coherence,[0,round(wd_shift*zp_fact_range_cnf)]);
end
end


