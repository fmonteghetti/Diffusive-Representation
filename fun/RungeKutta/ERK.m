function [q,t,Nop,T] = ERK(L,q_0,tf,dt,A,b)
%ERK Solve the system of ODEs:
%       q' = L(q,t) (dim. N, timestep dt)
% with initial condition on q (t=0).
% The system is solved through an explicit M-stage Runge-Kutta (ERK) method.
% Inputs:
%   L           function handle (Nx1,1x1) -> (Nx1)
%               ex: L(q,t) = A*q + u(t) (linear system with input u)
%   q_0 (Nx1) initial condition on q (t=0)
%   tf (1x1     final time
%   dt (1x1)   time-step
%   A (MxM)    Lower triangular matrix with null diagonal
%              (In Butcher notations: A(i,j) = a i,j)
%   B (M)      Vector of quadrature weights (In Butcher notations: bi)
% Outputs:
%   q1 (N1x(Nit+1)) computed values. q(:,i) is an approx. of q(t(i)).
%   t (1x(Nit+1)) corresponding instants : t(1) is 0, t(end) is tf.
%   Nop (1x1) Estimate of the number of operations performed.
%   T (1x1) Elapsed time
%
% Copyright (c) 2018 Florian Monteghetti.
% This work is licensed under the terms of the MIT license.  
% For a copy, see <https://opensource.org/licenses/MIT>.



    % Check validity of the arrays A and b.
    % A: lower triangular, null diagonal
validateattributes(A,{'numeric'},{'2d','nonempty','real','finite','nonnan'});
validateattributes(b,{'numeric'},{'vector','nonempty','real','finite','nonnan','numel',length(A)});
b = b(:); % vector column
c = sum(A,2); % ci coefficients (Butcher notations)
M = size(A,1);
fprintf('[ERK]: %d-stage explicit Runge-Kutta.\n',M);

dt_f = rem(tf,dt); % Final time-step (for the last iteration)
N_it = floor(tf/dt); % No. of iter. to perform with <dt>
    if(N_it==0)
        dt = dt_f;
    else
        dt = [dt*ones(1,N_it),dt_f];
    end
N_it = N_it+1; % Total number of iterations to perform

    % Initialization
N = length(q_0);
q = zeros(N,N_it+1); % Output: q(n) = q(t_{n-1})
q(:,1)=q_0;
t = zeros(1,N_it+1);
t(2:(end-1)) = dt(1:(end-1)).*(1:(N_it-1)); t(end)=t(end-1)+dt_f;
    % Estimate of the number of operations
Nop = N_it*N;
fprintf('[ERK]: %d iter. & %d oper. to perform (dt=%1.2g).\n',N_it,Nop,dt(1));
reverseStr='';
time = cputime();
    % Temporary values
k=zeros(N,M);
for n=1:N_it
   msg = sprintf('[ERK]: Percent done: %3.1f', floor(100*n/N_it));
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
   for i=1:M % M stage
       k(:,i) = L(q(:,n)+dt(n)*k*(A(i,:)'),t(n)+c(i)*dt(n));       
    end
    q(:,n+1)=q(:,n) + dt(n)*(k*b);
end
    T = cputime()-time;
    fprintf(' - Finished. (%2.2g s)\n',T);
end
