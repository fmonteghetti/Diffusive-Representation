function [A,b] = RK_convertCoeffs(A2N,B2N)
%RK_convertCoeffs Convert coefficient of a 2N-storage RK into the original 
% (A,b) (Butcher notations). Up to seven stages.
% Based on Appendix of 10.1016/j.jcp.2009.02.015.
% Rmk: accept symbolic or numerical inputs.
% Inputs:
%   A2N (N) Runge-Kutta coeffs, 2N-storage format.
%   B2N (N) Runge-Kutta coeffs, 2N-storage format.
% Outputs:
%   A (NxN) Runge-Kutta coeffs.
%   b (N)   Runge-Kutta quadrature weights.
%
% Copyright (c) 2018 Florian Monteghetti.
% This work is licensed under the terms of the MIT license.  
% For a copy, see <https://opensource.org/licenses/MIT>.


if length(A2N)>8
   error('[%s] Implemented formulas cover only up to 8 stages.',mfilename); 
end

N = length(A2N); % N-stage Explicit Runge-Kutta method
A = zeros(N,N);
b = zeros(N,1);

if N<8 % zero padding
   A2N((end+1):8) = 0;
   B2N((end+1):8) = 0;
end
    Ai = A2N;
    Bi = B2N;
    % Conversion to the original coefficients
    % Appendix of 10.1016/j.jcp.2009.02.015 provides formula for N<=7.
    % They have been extended to N=8.
        % quadrature weights
    b_tmp = zeros(8,1);
    b_tmp = [Ai(2)*(Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*(Ai(8)*Bi(8)+Bi(7))+Bi(6))+Bi(5))+Bi(4))+Bi(3))+Bi(2))+Bi(1);
            Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*(Ai(8)*Bi(8)+Bi(7))+Bi(6))+Bi(5))+Bi(4))+Bi(3))+Bi(2);
            Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*(Ai(8)*Bi(8)+Bi(7))+Bi(6))+Bi(5))+Bi(4))+Bi(3);
            Ai(5)*(Ai(6)*(Ai(7)*(Ai(8)*Bi(8)+Bi(7))+Bi(6))+Bi(5))+Bi(4);
            Ai(6)*(Ai(7)*(Ai(8)*Bi(8)+Bi(7))+Bi(6))+Bi(5);
            Ai(7)*(Ai(8)*Bi(8)+Bi(7))+Bi(6);
            (Ai(8)*Bi(8)+Bi(7));
            Bi(8)];
    for i=1:N
        b(i) = b_tmp(i);
    end
        % A matrix (lower triangular with null diagonal)
    A_tmp = zeros(8,8);
    A_tmp(2,1) =    [Bi(1)]; % OK
    A_tmp(3,1:2) =  [Ai(2)*Bi(2)+Bi(1),Bi(2)]; % OK
    A_tmp(4,1:3) =  [Ai(2)*(Ai(3)*Bi(3)+Bi(2))+Bi(1),Ai(3)*Bi(3)+Bi(2),Bi(3)]; % OK
    A_tmp(5,1:4) =  [Ai(2)*(Ai(3)*(Ai(4)*Bi(4)+Bi(3))+Bi(2))+Bi(1),Ai(3)*(Ai(4)*Bi(4)+Bi(3))+Bi(2),Ai(4)*Bi(4)+Bi(3),Bi(4)]; % OK
    A_tmp(6,1:5) =  [Ai(2)*(Ai(3)*(Ai(4)*(Ai(5)*Bi(5)+Bi(4))+Bi(3))+Bi(2))+Bi(1),Ai(3)*(Ai(4)*(Ai(5)*Bi(5)+Bi(4))+Bi(3))+Bi(2),Ai(4)*(Ai(5)*Bi(5)+Bi(4))+Bi(3),Ai(5)*Bi(5)+Bi(4),Bi(5)]; % OK
    A_tmp(7,1:6) =  [Ai(2)*(Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*Bi(6)+Bi(5))+Bi(4))+Bi(3))+Bi(2))+Bi(1),Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*Bi(6)+Bi(5))+Bi(4))+Bi(3))+Bi(2),Ai(4)*(Ai(5)*(Ai(6)*Bi(6)+Bi(5))+Bi(4))+Bi(3),Ai(5)*(Ai(6)*Bi(6)+Bi(5))+Bi(4),Ai(6)*Bi(6)+Bi(5),Bi(6)]; % OK    
    A_tmp(8,1:7) =  [Ai(2)*(Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5))+Bi(4))+Bi(3))+Bi(2))+Bi(1),Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5))+Bi(4))+Bi(3))+Bi(2),Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5))+Bi(4))+Bi(3),Ai(5)*(Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5))+Bi(4),Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5),Ai(7)*Bi(7)+Bi(6),Bi(7)]; % OK    
    for i=1:N
        A(i,:) = A_tmp(i,1:N);
    end
end
% Correct for N<=7
% b_tmp = [Ai(2)*(Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5))+Bi(4))+Bi(3))+Bi(2))+Bi(1);
%         Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5))+Bi(4))+Bi(3))+Bi(2);
%         Ai(4)*(Ai(5)*(Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5))+Bi(4))+Bi(3);
%         Ai(5)*(Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5))+Bi(4);
%         Ai(6)*(Ai(7)*Bi(7)+Bi(6))+Bi(5);
%         Ai(7)*Bi(7)+Bi(6);
%         Bi(7)];
% A_tmp(2,1) =    [Bi(1)]; % OK
% A_tmp(3,1:2) =  [Ai(2)*Bi(2)+Bi(1),Bi(2)]; % OK
% A_tmp(4,1:3) =  [Ai(2)*(Ai(3)*Bi(3)+Bi(2))+Bi(1),Ai(3)*Bi(3)+Bi(2),Bi(3)]; % OK
% A_tmp(5,1:4) =  [Ai(2)*(Ai(3)*(Ai(4)*Bi(4)+Bi(3))+Bi(2))+Bi(1),Ai(3)*(Ai(4)*Bi(4)+Bi(3))+Bi(2),Ai(4)*Bi(4)+Bi(3),Bi(4)]; % OK
% A_tmp(6,1:5) =  [Ai(2)*(Ai(3)*(Ai(4)*(Ai(5)*Bi(5)+Bi(4))+Bi(3))+Bi(2))+Bi(1),Ai(3)*(Ai(4)*(Ai(5)*Bi(5)+Bi(4))+Bi(3))+Bi(2),Ai(4)*(Ai(5)*Bi(5)+Bi(4))+Bi(3),Ai(5)*Bi(5)+Bi(4),Bi(5)]; % OK
% A_tmp(7,1:6) =  [Ai(2)*(Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*Bi(6)+Bi(5))+Bi(4))+Bi(3))+Bi(2))+Bi(1),Ai(3)*(Ai(4)*(Ai(5)*(Ai(6)*Bi(6)+Bi(5))+Bi(4))+Bi(3))+Bi(2),Ai(4)*(Ai(5)*(Ai(6)*Bi(6)+Bi(5))+Bi(4))+Bi(3),Ai(5)*(Ai(6)*Bi(6)+Bi(5))+Bi(4),Ai(6)*Bi(6)+Bi(5),Bi(6)]; % OK    
