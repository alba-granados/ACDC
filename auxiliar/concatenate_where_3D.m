%% HEADER
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% GPP 
% This code concatenates to vector given the order where you want them.
%
% 
% ----------------------------------------------------------
% Author:    Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/09/2013)
%
% Version  record
% 1.0 2016/03/17    Imported code from S6

function [X]= concatenate_where_3D (x,y,where)

    if(strcmp(where,'back'))
        X = [y ; x];
    elseif(strcmp(where,'fore'))
        X = [x ; y];
    end
    
end