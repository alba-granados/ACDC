% calculate the window size for each stack in order to fit all the beams
% without any wrapp

function [L1BS]= avoid_wrapps(L1BS)

global N_samples_sar_chd N_ku_pulses_burst_chd

L1BS.N_windows = ceil(max(max(L1BS.shift)/N_samples_sar_chd));
truncated_zeros    = zeros(max(L1BS.N_beams_stack),N_samples_sar_chd*(L1BS.N_windows-1));
for i_surf = 1:L1BS.N_total_surf_loc
    beams_surf_aux1_tmp             = squeeze(L1BS.beams_surf(i_surf,:,:));
    beams_surf_aux1_fft             = fftshift(fft(beams_surf_aux1_tmp.'),1).';
    beams_surf_aux1_fft_increased   = (cat(2,beams_surf_aux1_fft,truncated_zeros));
    
    beams_surf_aux1_tmp_increased   = ifft(ifftshift(beams_surf_aux1_fft_increased(:,1:L1BS.N_windows:end),2).').';

    L1BS.beams_surf(i_surf,:,:)     = beams_surf_aux1_tmp_increased;
    
    %% es perd resoluci'o al fer al ifft amb menys mostres
   
    
    
    
    
end

end