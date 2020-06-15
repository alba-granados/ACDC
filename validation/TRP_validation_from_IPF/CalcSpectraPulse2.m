% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
% 
% This code implements the algorithm as described in the
% ISARD_ESA_CR2_TRP_CAL_DPM_030 2.b of 26/05/2011
%
% ---------------------------------------------------------
%
% CALCSPECTRAPULSE2: Calculates the fourier transform. Zero padding is optional
% 
% Calling
%   Sy = CalcSpectraPulse2 (y, paddingfactor)
%
% INPUTS
%============
%   y : input signal. Column vector !!
%   paddingfactor: a padding factor of, for example, 64 means that 63 zeros are inserted every gate.
%
% OUTPUT
%============
%  Sy: fourier transform of y, complex number.Column vector.
%
% Comments
%   This program works for only one echoe.
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
%           Daniel Martinez / Pildo Labs
%           Josep Montolio / Pildo Labs
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Josep Montolio / Pildo Labs (20/10/09)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Sy=CalcSpectraPulse2(y, paddingfactor)

y=flippad(y, paddingfactor);

N=length(y);
Sy=1/N.*fftshift(fft(y)); % Normalized fft

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol=flippad(y, paddingfactor);
% Flips and pad with zeros
% A padding factor of 64 means that 63 zeros are inserted every gate.

L=length(y);
flip = y;

sol=[   flip(1:(L/2),1)
        zeros((paddingfactor-1)*L,1)
        flip((L/2+1):L,1)];
 


