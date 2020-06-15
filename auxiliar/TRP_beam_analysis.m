function [central, central_beams] = TRP_beam_analysis(L1BS)
global pi_cst N_samples 
zp_fact_trp = 1;
    central_beams = 0;
    beam_power = zeros (1,L1BS.N_beams_stack);
    for i_beam = 1:L1BS.N_beams_stack
        beam_power(i_beam) = sum (abs(fftshift(fft(L1BS.beam_geo_corr(i_beam,:).',N_samples*zp_fact_trp),1)).'.^2); 
    end
    beam_power=beam_power/max(beam_power);
    [~,pos_max]=max(beam_power);
    
%     [~,doppler_central]= min(abs(pi_cst/2-L1BS.beam_ang_surf(:)));
    right_beams_cnf=50;
    if(pos_max<50)
        left_beams_cnf=pos_max-1;
    else
        left_beams_cnf=50;
    end

    if((pos_max+right_beams_cnf)>L1BS.N_beams_stack)
        right_beams_cnf=L1BS.N_beams_stack-pos_max;
    end
    central_beams_index = (pos_max-left_beams_cnf:pos_max+right_beams_cnf);
    valid_beams     = find(L1BS.Gap_flag(pos_max-left_beams_cnf:pos_max+right_beams_cnf));
    central_beams = central_beams_index(valid_beams);
    central_beams = central_beams(beam_power(central_beams)>0.1);
    try
    % Gaussian fittting
         cfun_beams                       = fit ((central_beams)', beam_power(central_beams)', 'gauss1'); %fitting vs look_angle
%         
        a_gauss                 = cfun_beams.a1;
        stack_beam_centre       = cfun_beams.b1;
        stack_std               = cfun_beams.c1/2;
        stack_width             = 2*sqrt(2*log(2))*cfun_beams.c1;
        % Compute the characterization parameters:
%         stack_skewness= skewness(power_fitted(central_beams));
%         stack_kurtosis= kurtosis(power_fitted(central_beams))-3;

        cfun_beams = fit ((central_beams)', beam_power(central_beams)', 'gauss1'); %fitting vs look_angle
        stack_beam_centre       = cfun_beams.b1;
    catch
        stack_look_ang_centre       = 0;
        stack_pointing_ang_centre   = 0;
        stack_std               = 0;
        stack_width             = 0;
        stack_skewness= 0;
        stack_kurtosis= 0;
        stack_beam_centre       = 0;

    end
central =  (round(stack_beam_centre));
end
