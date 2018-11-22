function [A2N,B2N,C] = RK_get84Coeffs_2NStorage()
%RK_get84Coeffs_2NStorage 
%Return the coefficients of a 8-stage 4th order Runge-Kutta (RKF84)
%method suited for aeroacoustics with D.G. (10.1016/j.jcp.2011.11.024).
% The scheme is given in 2N-storage format, with 16-digit precision.
% Outputs:
%   A2N (5)     "A" coefficients
%   B2N (5)     "B" coefficients
%   C (1)       Stability region: |lambda*dt|<C
%
% Copyright (c) 2018 Florian Monteghetti.
% This work is licensed under the terms of the MIT license.  
% For a copy, see <https://opensource.org/licenses/MIT>.

    warning('[%s] Coefficients only known with 16-digit precision',mfilename);
    A2N = [0,-0.5534431294501569,0.01065987570203490,-0.5515812888932000,-1.885790377558741,-5.701295742793264,2.113903965664793,-0.5339578826675280];
    B2N = [0.08037936882736950,0.5388497458569843,0.01974974409031960,0.09911841297339970,0.7466920411064123,1.679584245618894,0.2433728067008188,0.1422730459001373];
    C = 7.9;
end
