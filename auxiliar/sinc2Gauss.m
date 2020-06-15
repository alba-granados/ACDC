function [A_s2G,alpha_g] = sinc2Gauss(alpha_cnf, n)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% This code generates the configuration file input to isardSATs
% re-tracker as described in isardSAT_Retracker_DPM_V1.a.docx internal
% document
% 
%
% ---------------------------------------------------------
% sinc2Gauss: This code measures the sinc to Gaussian approximation in terms of alpha_g
% parameter. Note that
%
% sinc?2(m) ~ exp^(alpha_g*m?2)
%
% Calling
%   sinc2Gauss
%
% Inputs
%   alpha_cnf   -   Apodization function alpha parameter
%   n - rank of values it determines number of samples
%   
%   window (n) = alpha_cnf + (1-alpha_cnf) cos(2 pi n / (N-1)) n = 1...N,
%   or 0 ... N-1

%   if alpha_cnf = 1 --> boxcard
%   if alpha_cnf = 0.5 --> Hanning
%   if alpha_cnf = 0.54 --> Hamming
%
% Output
%   alpha_g parameter to scale Gaussian fon
%
%
% ----------------------------------------------------------
% 
% Author:   Cristina Martin-Puig / isardSAT
%
% Reviewer: Cristina Martin-Puig / isardSAT
%
% Last revision: Cristina Martin-Puig / isardSAT 25/02/2013
%
%
% This software is subject to the conditions 
% set forth in the ESA contract P4 Jason-CS GPP
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define initial values for fitting parameters
x = [0.5 0.5]; % alpha_g and A gaussian initial values
%x=0.5;
% Calculate Dataset to be fitted
window      =   alpha_cnf*DeltaN(n)-(1-alpha_cnf)/2*(DeltaN(n-length(n)/(length(n)-1))+DeltaN(n+length(n)/(length(n)-1)));
window      =   abs(window).^2;

% set fitting option
options     =   optimset('Algorithm', 'levenberg-marquardt','Display','off');


% fit real with theoretical 
x =   lsqcurvefit (@Gaussian,x,n,window,[ ],[ ],options);

% Assigne results to data 
A_s2G       =   x(1); % Amplitude from sinc to Gaussian conversion
alpha_g     =   x(2); % alpha_g parameter

% plot (n,window,'b');
% hold on
% plot (n,Gaussian(x,n),'m');
% grid on
% title (['MLE result; alpha-g = ',num2str(alpha_g),' Amp = ', num2str(A_s2G)])
% ylabel (['MLE fitting of along-track SPTR with N = ',num2str(length(n)),' samples'])
% xlabel ('Doppler bin resolution units')

% 
% print('-dpdf','isardSAT_SARMret_alphag.pdf')
% close all;

end




