function IRF_sinc2gaussian(N_samples,type_window)


resampling_factor=16; %useful for well approximation

IRF_sinc_ref=fft(ones(1,N_samples),resampling_factor*N_samples);
x_vec=-resampling_factor*N_samples/2:(resampling_factor*N_samples/2-1);
switch type_window
    case 'boxcar'
        window=ones(1,N_samples);
        
    case 'hamming'
        window=hamming(N_samples);
    case 'hanning'
        window=hann(N_samples);
        
end
IRF_wind=(fft(window,resampling_factor*N_samples)/(sqrt(resampling_factor*N_samples)*mean_window)).^2;
init_param=[]

% set fitting option
options     =   optimset('Algorithm', 'levenberg-marquardt','Display','off');

    function res=gaussian(x_vec,A_g,sigma_g)
        res=A_g*exp(-x_vec.^2/(2.0*sigma_g^2));
    end
% fit real with theoretical 
x =   lsqcurvefit (@gaussian,,x_vec,IRF_wind,[],[],options);


end