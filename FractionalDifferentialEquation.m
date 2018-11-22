%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solution of the scalar fractional differential equation:
%            dy/dt(t) = a*y(t) - g * (d^alpha)y(t) + u(t) (t>0)
%            y(0) = y0,
% where d^alpha is the Caputo derivative.
%
% The numerical solution is based on a discretization of the diffusive
% representation of the fractional derivative through a quadrature-based 
% method, eventually refined with a nonlinear least squares optimization.
%
% Outline of the script:
%   I   - Compute discrete diffusive representation(s)
%   II  - Time-integration
%   III - Analytical solution for alpha=0.5 and u(t)=0
%   IV  - Plot solution(s)
%   V   - Export to CSV file
%
% Script tested on MATLAB R2015b
% Copyright (c) 2018 Florian Monteghetti.
% This work is licensed under the terms of the MIT license.  
% For a copy, see <https://opensource.org/licenses/MIT>.
    % FDE parameters
y0=1;
a = -1;
g = 1*exp(1i*(0*pi));
u=@(t)(0*cos(t));
    % Cell storing discrete diffusive representation
FracDerivative = cell(0);
N = 6; % Number of quadrature nodes
alpha = 0.5; % order of fractional derivative
mu_an = @(xi)(sin(alpha*pi)./(pi*xi.^alpha)); % diffusive weight (xi>0)
beta = @(alpha)min(alpha,1-alpha); % change of variable parameter
%% I - Build diffusive representation (quadrature rule alone)
[xi_quad,mu_quad]=ComputeDiscreteDiffusiveRep_Quadrature(mu_an,N,'CoVParam',beta(alpha));
FracDerivative{end+1}=struct('xi',xi_quad,'mu',mu_quad,'Name',sprintf('Quad-min-ximin=%.3e-ximax=%.3e-Nxi=%d',min(xi_quad),max(xi_quad),N));
%% I - Build diffusive representation (quadrature then nonlinear optimization)
% Standard cost function
Tol = 1e-15; % Termination tolerance
K=1e4; % Number of angular frequencies in cost function
[xi_quad,mu_quad]=ComputeDiscreteDiffusiveRep_Quadrature(mu_an,N,'CoVParam',beta(alpha),'Optim',struct('Tol',Tol,'K',K,'Method','NonLinOptimStd','h_an',@(s)s.^(-alpha)));
FracDerivative{end+1}=struct('xi',xi_quad,'mu',mu_quad,'Name',sprintf('Quad-Min-Opt-Std-Log-NL-Tol=%1.3e-K=%d-ximin=%.3e-ximax=%.3e-N=%d',Tol,K,min(xi_quad),max(xi_quad),N));
%% I - Build diffusive representation (quadrature then nonlinear optimization)
% Standard cost function, normalized
Tol = 1e-15; K=1e4;
[xi_quad,mu_quad]=ComputeDiscreteDiffusiveRep_Quadrature(mu_an,N,'CoVParam',beta(alpha),'Optim',struct('Tol',Tol,'K',K,'Method','NonLinOptimStdNorm','h_an',@(s)s.^(-alpha)));
FracDerivative{end+1}=struct('xi',xi_quad,'mu',mu_quad,'Name',sprintf('Quad-Min-Opt-StdNorm-Log-NL-Tol=%1.3e-K=%d-ximin=%.3e-ximax=%.3e-N=%d',Tol,K,min(xi_quad),max(xi_quad),N));
%% I - Build diffusive representation (quadrature then nonlinear optimization)
% Extended cost function
Tol = 1e-15; K=1e4;
[xi_quad,mu_quad]=ComputeDiscreteDiffusiveRep_Quadrature(mu_an,N,'CoVParam',beta(alpha),'Optim',struct('Tol',Tol,'K',K,'Method','NonLinOptimExt','h_an',@(s)s.^(-alpha)));
FracDerivative{end+1}=struct('xi',xi_quad,'mu',mu_quad,'Name',sprintf('Quad-Min-Opt-Ext-Log-NL-Tol=%1.3e-K=%d-ximin=%.3e-ximax=%.3e-N=%d',Tol,K,min(xi_quad),max(xi_quad),N));
%% I - Build diffusive representation (quadrature then nonlinear optimization)
% Reflection coefficient cost function
Tol = 1e-15; K=1e4;
[xi_quad,mu_quad]=ComputeDiscreteDiffusiveRep_Quadrature(mu_an,N,'CoVParam',beta(alpha),'Optim',struct('Tol',Tol,'K',K,'Method','NonLinOptimRefl','h_an',@(s)s.^(-alpha)));
FracDerivative{end+1}=struct('xi',xi_quad,'mu',mu_quad,'Name',sprintf('Quad-Min-Opt-Refl-Log-NL-Tol=%1.3e-K=%d-ximin=%.3e-ximax=%.3e-N=%d',Tol,K,min(xi_quad),max(xi_quad),N));
%% I - Build diffusive representation (quadrature then nonlinear optimization) 
% Reflection coefficient cost function, normalized
Tol = 1e-15; K=1e4;
[xi_quad,mu_quad]=ComputeDiscreteDiffusiveRep_Quadrature(mu_an,N,'CoVParam',beta(alpha),'Optim',struct('Tol',Tol,'K',K,'Method','NonLinOptimReflNorm','h_an',@(s)s.^(-alpha)));
xi_min  = min(xi_quad); xi_max = max(xi_quad);
FracDerivative{end+1}=struct('xi',xi_quad,'mu',mu_quad,'Name',sprintf('Quad-Min-Opt-ReflNorm-Log-NL-Tol=%1.3e-K=%d-ximin=%.3e-ximax=%.3e-N=%d',Tol,K,min(xi_quad),max(xi_quad),N));
%% II - Time integration with Explicit Runge Kutta (ERK)
y_sol = cell(0); % cell storing computed solution(s)
dt_ERK = 9e-03; % time step
tf=100;
for i=1:length(FracDerivative) % for each discrete representation
    N = length(FracDerivative{i}.xi);
        % Build coupled system
    A = zeros(1+N);
    A(1,1) = a-g*sum(FracDerivative{i}.mu(:)); A(1,2:end)=g*transpose(FracDerivative{i}.mu(:).*FracDerivative{i}.xi(:));
    A(2:end,1) = 1; A(2:end,2:end)=-diag(FracDerivative{i}.xi(:));
        % Get RK coeffs
    [A_ERK,b_ERK] = RK_get84Coeffs_2NStorage();
    [A_ERK,b_ERK] = RK_convertCoeffs(A_ERK,b_ERK);
        % Initial condition
    Y0 = y0*[1;1./FracDerivative{i}.xi(:)];
        % Time integration
    [y_diff,t_diff] = ERK(@(y,t)(A*y+[u(t);zeros(N,1)]),Y0,tf,dt_ERK,A_ERK,b_ERK);
        % Extract solution
    y_diff =  y_diff(1,:);
    y_sol{i} = [t_diff(:),y_diff(:)];
end
%% III - Analytical solution of the FDE for u(t)=0 and alpha=0.5
% Solution of the scalar fractional differential equation:
%           (d^1/2)^2 y(t) = a*y(t) - g * (d1/2) y(t) + u(t) (t>0)
%           y(0) = y0, y'(0) = y1,
% where d1/2 is the Caputo derivative.
% To compute the analytical solution, this is rewritten as
% P(d1/2)y=u, with P(s)=(s-l1)(s-l2):
%           (d1/2)^2 y(t) -(l1+l2)d1/2 y(t) + l1*l2*y(t) = u(t) (t>0)
% See D. Matignon, An introduction to fractional calculus, in: Scaling, 
% Fractals and Wavelets, ISTE–Wiley, London–Hoboken, 
% ISBN 978-1-84-821072-1, pp. 237–277, doi:10.1002/9780470611562.ch7, 2009.
y1=0; % If y1=0,  then (d1/2)^2 = y'(t).
        % Analytical solution
        % Mittag-Leffler functions E_(1/2),1
l1 = roots([1,g,-a]); l1=l1(1); l2 = -g-l1;
E_ex = @(l,t)(exp((l^2)*t).*(1+erf_(l*sqrt(t)))); % ((1x1),(1xm) -> (1xm))
y_ex = @(t)( (y0/(l1-l2))*(l1*E_ex(l2,t)-l2*E_ex(l1,t)) + (y1/(l1-l2))*(E_ex(l1,t)-E_ex(l2,t)));
        % Equivalent alternative using code by 
        % Roberto Garrappa, University of Bari
        % https://fr.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function
%E_ex = @(l,t)MittagLeffler(l*sqrt(t),0.5,1);
%y_ex = @(t)( (y0/(l1-l2))*(l1*E_ex2(l2,t)-l2*E_ex2(l1,t)) + (y1/(l1-l2))*(E_ex2(l1,t)-E_ex2(l2,t)));
clear l1 l2
%% IV - Plot computed solution
t=linspace(0,tf,1e2);
figure
clf
leg=cell(0);
subplot(2,1,1)
hold all
for i=1:length(y_sol) % plot computed solution(s)
    plot(y_sol{i}(:,1),y_sol{i}(:,2));
    leg(end+1)={FracDerivative{i}.Name};
end
if alpha==0.5 % plot analytical solution
    plot(t,y_ex(t),'k--');
    leg(end+1)={sprintf('Exact')};
end
legend(leg);
title('Fractionnal differential equation');
xlabel('t');
ylabel('y');
ylim([0,1])
subplot(2,1,2)
hold all
for i=1:length(y_sol)
    t_err = y_sol{i}(:,1);
    plot(t_err,100*(y_ex(t_err)-y_sol{i}(:,2))./y_ex(t_err),'-');
end
ylabel('Relative error (%)');
xlabel('t');
ylim([-10,10])
%% V - Export exact solution
% CSV file with columns:
%   t,y
t = linspace(0,tf,1e2);
csvhead = sprintf('t,y');
namestr = sprintf('FDE_a=%.3e_g=%.3e_y0=%.3e_tf=%.3e_ExactSol.csv',a,g,y0,tf);
dlmwrite(namestr,csvhead,'Delimiter','');
dlmwrite(namestr,[t(:),y_ex(t)'],'-append','Delimiter',',','newline','unix','precision','%1.6e');
%% V - Export computed solution
% CSV file with columns:
%   t,y
for i=1:length(y_sol)
    t = linspace(0,tf,1e2);
    csvhead = sprintf('t,y');
    namestr = sprintf('FDE_a=%.3e_g=%.3e_y0=%.3e_tf=%.3e_CompSol_dt=%.3e_%s.csv',a,g,y0,tf,dt_ERK,FracDerivative{i}.Name);
    dlmwrite(namestr,csvhead,'Delimiter','');
    t_csv = y_sol{i}(:,1); 
    idx=1:ceil(length(t_csv)/300):length(t_csv);
    t_csv = t_csv(idx);
    y_csv = y_sol{i}(idx,2);
    dlmwrite(namestr,[t_csv(:),y_csv],'-append','Delimiter',',','newline','unix','precision','%1.6e');
end
%% V - Export relative error
% a,g,y0, tf, dt
csvhead = sprintf('t,RelError(percent)');
for i=1:length(y_sol)
    namestr = sprintf('FDE_a=%.3e_g=%.3e_y0=%.3e_tf=%.3e_RelError_dt=%.3e_%s.csv',a,g,y0,tf,dt_ERK,FracDerivative{i}.Name);
    dlmwrite(namestr,csvhead,'Delimiter','');
    t_err = y_sol{i}(:,1);
    y_err = 100*(y_ex(t_err)-y_sol{i}(:,2))./y_ex(t_err);
        % undersample (length around 100)
    idx=1:ceil(length(t_err)/300):length(t_err);
    y_err = y_err(idx); t_err=t_err(idx);
    dlmwrite(namestr,[t_err,y_err],'-append','Delimiter',',','newline','unix','precision','%1.6e');
end
