%% save different images for the same plot with different angles "view"


%  rotateFigure(h,axis)
lin=1;
log=0;
velocity='ecc';
 N_images=250; % N_images/25 = video duration with reasonable quality
 elev_start=0;
 elev_end=90;
 azimut_start=0;
 azimut_end = 180;
 switch velocity
    case 'lin'
        az = linspace(azimut_start,azimut_end,N_images);
        el = linspace(elev_start,elev_end,N_images);
     case 'log'
        az = logspace(log10(azimut_end),log10(azimut_start),N_images);
        el = logspace(log10(elev_start),log10(elev_end),N_images);
     case 'ecc'
        velocity_factor=1; % values between 0.5 to 2 are fine
        az = (((linspace(-azimut_end,azimut_end,N_images)).^(2^velocity_factor))./(azimut_end).^(2^velocity_factor-1));
        az(1:N_images/2)=-az(1:N_images/2);
        az=az/2+azimut_end/2;
        el = (((linspace(-elev_end,elev_end,N_images)).^(2^velocity_factor))./(elev_end).^(2^velocity_factor-1));
        el(1:N_images/2)=-el(1:N_images/2);
        el=el/2+elev_end/2;
 end    
 for i_view= 1:N_images
      view(az(i_view),el(i_view));
    
%      pause(1);
     saveas(gcf, sprintf('%03i_sample.png',(i_view)) , 'png');
 end