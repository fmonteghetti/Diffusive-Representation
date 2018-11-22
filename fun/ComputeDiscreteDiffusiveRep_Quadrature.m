function [xi,mu] = ComputeDiscreteDiffusiveRep_Quadrature(mu_an,N,varargin)
%ComputeDiscreteDiffusiveRep_Quadrature Compute a discretisation of the
% diffusive representation [xi,mu]:
%      -
%     / mu_an(xi)                             muk   
%    -  --------- dxi (continuous) -> Sum   ------- (discrete)
%        s + xi                             s + xik
% The discrete poles and weights are computed using a Gauss-Legendre 
% quadrature rule on the interval (-1,1) combined with a change of variable
% parametrized by the nonnegative scalar beta.
% Optionally, the computed weights [xi,mu] can be refined using a nonlinear
% least squares optimization with constraints (MATLAB function lsqcurvefit): 
% several cost functions are available. The constraints are:
%               muk>=0, xik>=0, xik<=xi_max.
% Created on MATLAB R2015b.
%  
% Note: The change of variable and cost functions are described in 
%   F. Monteghetti, D. Matignon, E. Piot. Time-local discretization of 
%   fractional and related diffusive operators using Gaussian quadrature 
%   with applications (Applied Numerical Mathematics, 2018)
%
% Inputs:
%   N (1) number of quadrature nodes
%   mu_an (function handle) analytical expression of the diffusive weight
% Inputs (optional):
%   CoVParam (1) "beta" parameter for the change of variable (>=0)
%                (default: 1/2)
%   Optim (struct) parameter that describes the optimization to be used to 
%                  refine poles and weights.
%                  Fields:
%       - h_an: transfer function @(s)h_an(s) to be approximated
%       - Method: optimization method
%           'NonLinOptimStd': nonlinear least squares with cost
%           function J.
%           'NonLinOptimStdNorm': nonlinear least squares with cost
%           function J, normalized.
%           'NonLinOptimExt': nonlinear least squares with cost
%           function Js.
%           'NonLinOptimRefl': nonlinear least squares with cost
%           function Jpsi.
%           'NonLinOptimReflNorm': nonlinear least squares with cost
%           function Jpsi, normalized.
%       - (Optional) K: number of angular frequencies used to define cost
%                    function. (default: 1e2)
%       - (Optional) Tol: tolerance (x & fun) for optimization
%                    (default: 1e-15)
% Outputs:
%   xik (N) poles
%   muk (N) weights
%
% Copyright (c) 2018 Florian Monteghetti.
% This work is licensed under the terms of the MIT license.  
% For a copy, see <https://opensource.org/licenses/MIT>.


        % Validate input arguments
    validateattributes(mu_an,{'function_handle'},{},mfilename,'mu_an');
    xi = [0,1];
    validateattributes(mu_an(xi),{'numeric'},{'size',size(xi)},mfilename,'mu_an(xi)');
    validateattributes(N,{'numeric'},{'scalar','integer','positive'},mfilename,'N');
        % Retrieve and validate optional arguments
        % ! May not work on older MATLAB versions. (If so, remove the
        % verifications.)
    p = inputParser; p.FunctionName = mfilename;
    addParameter(p,'CoVParam',0.5,@(x)(validateattributes(x,{'numeric'},{'scalar','positive'})));
    addParameter(p,'Optim',[],@(x)validateOptimStruct(x));
    parse(p,varargin{:});  oA = p.Results;
    beta = oA.CoVParam; % parameter for change of variable
    Optim=oA.Optim; % structure describing optimization
            % Get quadrature nodes on (-1,1) (no end points)
    [x,w] = GaussLegendreQuadraturePoints(N);
            % Compute diffusive weights and poles
    xi = (1+x)./(1-x); xi = xi.^(1/beta);
    mu = (1-x).^(-1-1/beta).*(1+x).^(-1+1/beta).*mu_an(xi);
    mu = 2/beta*mu.*w;
    fprintf('[%s] Computation of weights/poles with quadrature (beta=%1.3e)....\n',mfilename,beta);
    fprintf('[%s] Done. xi_max=%1.3e, xi_min=%1.3e, mu_max=%1.3e, mu_in=%1.3e.\n',mfilename,max(xi),min(xi),max(mu),min(mu));
        % Refinment with nonlinear least squares optimization
    if isempty(Optim) 
        return; % quit now if no optimization is requested
    end
    
    h_an = Optim.h_an;
        % Number of angular frequencies in cost function
    K = 1e2*length(xi);
    if isfield(Optim,'K')
       K = Optim.K;
    end
    omegan = logspace(log10(min(xi)),log10(max(xi)),K);
        % Termination tolerance for optimization
    Tol=1e-15;
    if isfield(Optim,'Tol')
       Tol = Optim.Tol;
    end
    
    N = length(xi);
        % Definition of coefficients for optimization
    coeff_init = [xi(:);mu(:)];
    coeff2xi = @(coeff)coeff(1:N);
    coeff2mu = @(coeff)coeff((N+1):end);
    options = optimoptions('lsqcurvefit','Jacobian','on','DerivativeCheck','off','Display','final-detailed','TolX',Tol,'TolFun',Tol,'MaxIter',5e4','MaxFunEvals',5e4);
        % Lower bound constraints
            %   [xi,mu] >=0
    lb=0*coeff_init; % length 2*Nxi
        % Upper bound constraints
            %   xi<=xi_max, no constraints on mu
    ub = ones(N,1)*max(coeff2xi(coeff_init)); % length Nxi only
    
        % Choice of the cost function
    if strcmp(Optim.Method,'NonLinOptimStd')
       fprintf('[%s] Refinment with nonlinear least squares (Cost function J)...\n',mfilename);
       ydata = h_an(1i*omegan); ydata = [real(ydata(:));imag(ydata(:))];
       CostFun=@CostFunStandard;
   elseif strcmp(Optim.Method,'NonLinOptimStdNorm')
        fprintf('[%s] Refinment with nonlinear least squares (Cost function J, normalized)...\n',mfilename);
        ydata = [ones(1*length(omegan),1);zeros(1*length(omegan),1)];
        CostFun=@(coeff,xdata)CostFunStandardNormalized(coeff,xdata,h_an);
    elseif strcmp(Optim.Method,'NonLinOptimExt')
        fprintf('[%s] Refinment with nonlinear least squares (Cost function Js)...\n',mfilename);
        ydata = 1i*omegan(:).*h_an(1i*omegan(:)); ydata = [real(ydata(:));imag(ydata(:))];
        CostFun=@CostFunExtended;
   elseif strcmp(Optim.Method,'NonLinOptimRefl')
        fprintf('[%s] Refinment with nonlinear least squares (Cost function Jpsi)...\n',mfilename);
        ydata = z2beta(h_an(1i*omegan)); ydata = [real(ydata(:));imag(ydata(:))];
        CostFun=@CostFunRefl;
   elseif strcmp(Optim.Method,'NonLinOptimReflNorm')
        fprintf('[%s] Refinment with nonlinear least squares (Cost function Jpsi, normalized)...\n',mfilename);
        ydata = [ones(1*length(omegan),1);zeros(1*length(omegan),1)];
        CostFun=@(coeff,xdata)CostFunReflNormalized(coeff,xdata,h_an);
    else
        error('Unknown parameter for ''Optim.Method''.');
        return;
    end
        % Optimization
    [coeff_optim,~,~,exitflag] = lsqcurvefit(CostFun,coeff_init,omegan,ydata,lb,ub,options);
        % Retrieve xi and mu
    if exitflag<=0 % algorithm did not converge properly
        xi = NaN*coeff_optim;
        mu= NaN*coeff_optim;
    else
        xi = coeff2xi(coeff_optim);
        mu = coeff2mu(coeff_optim);
    end
end

function validateOptimStruct(x)
%validateOptimStruct Check wether x is a valid structure to describe the
%optimization stage.
% Input:
%   x (struct)
    if ~isstruct(x)
       error('Optim is not a structure.');
    end
    if ~isfield(x,'Method') || ~isfield(x,'h_an')
        error('Structure ''Optim'' has no field ''Method'' or ''h_an''.'); 
    end
end

% List of cost functions
function [F,J] = CostFunStandard(coeff,xdata)
%fun Function to be used in call to lsqcurvefit.
% Cost function J.
% Outputs:
%   F (2*length(xdata)) function values at coeff
%   J (2*length(xdata) x length(coeff)) Jacobian matrix at coeff

    K = length(xdata); % number of angular frequencies
    Nxi = length(coeff)/2; % number of poles
        % Function definition for optimization
    coeff2xi = @(coeff)coeff(1:Nxi);
    coeff2mu = @(coeff)coeff((Nxi+1):end);
    H_optim = @(s,coeff)get1stOrderUnitaryFilter(s,coeff2xi(coeff),'type','low-pass')*coeff2mu(coeff);
    F = [real(H_optim(1i*xdata(:),coeff));imag(H_optim(1i*xdata(:),coeff))];
if nargout > 1   % Two output arguments
        % Matrix: (1i*omega_i+xi_j)^-1 (size K X Nxi)
    H = get1stOrderUnitaryFilter(1i*xdata(:),coeff2xi(coeff),'type','low-pass');
        % Matrix: -mu_j/(1i*omega_i+xi_j)^2 (size K x Nxi)
    G = -(H.*H).*kron(ones(length(H),1),transpose(coeff2mu(coeff))); % a bit costly in memory
        % Jacobian matrix
    J = zeros(2*K,2*Nxi);
            % Real part w.r.t. xi (K*Nxi block)
    J(1:K,1:Nxi) = real(G);
            % Real part w.r.t. mu (K*Nxi block)
    J(1:K,Nxi+(1:Nxi)) = real(H);
            % Imag part w.r.t. xi (K*Nxi block)
    J(K+(1:K),1:Nxi) = imag(G);
            % Imag part w.r.t. mu (K*Nxi block)
    J(K+(1:K),Nxi+(1:Nxi)) = imag(H);
end
end

function [F,J] = CostFunStandardNormalized(coeff,xdata,h_an)
%fun Function to be used in call to lsqcurvefit.
% Cost function J, normalized.
% Outputs:
%   F (2*length(xdata)) function values at coeff
%   J (2*length(xdata) x length(coeff)) Jacobian matrix at coeff
%   h_an (fun_handle) Standard diffusive function being approximated (Laplace)

    K = length(xdata); % number of angular frequencies
    Nxi = length(coeff)/2; % number of poles
        % Function definition for optimization
    coeff2xi = @(coeff)coeff(1:Nxi);
    coeff2mu = @(coeff)coeff((Nxi+1):end);
    H_optim = @(s,coeff)get1stOrderUnitaryFilter(s,coeff2xi(coeff),'type','low-pass')*coeff2mu(coeff);
    tmp = H_optim(1i*xdata(:),coeff)./h_an(1i*xdata(:));
    F = [real(tmp);imag(tmp)];
if nargout > 1   % Two output arguments
        % Matrix: (1i*omega_i+xi_j)^-1 (size K X Nxi)
    H = get1stOrderUnitaryFilter(1i*xdata(:),coeff2xi(coeff),'type','low-pass');
        % Matrix: -mu_j/(1i*omega_i+xi_j)^2 (size K x Nxi)
    G = -(H.*H).*kron(ones(length(H),1),transpose(coeff2mu(coeff))); % a bit costly in memory
        % Mult matrix [1/H(1i*omega_i)]_i,j (size K x Nxi)
    Mult = kron(1./h_an(1i*xdata(:)),ones(1,Nxi));
        % Jacobian matrix
    J = zeros(2*K,2*Nxi);
            % Real part w.r.t. xi (K*Nxi block)
    J(1:K,1:Nxi) = real(Mult.*G);
            % Real part w.r.t. mu (K*Nxi block)
    J(1:K,Nxi+(1:Nxi)) = real(Mult.*H);
            % Imag part w.r.t. xi (K*Nxi block)
    J(K+(1:K),1:Nxi) = imag(Mult.*G);
            % Imag part w.r.t. mu (K*Nxi block)
    J(K+(1:K),Nxi+(1:Nxi)) = imag(Mult.*H);
end
end

function [F,J] = CostFunExtended(coeff,xdata)
%fun Function to be used in call to lsqcurvefit.
% Cost function Js.
% Outputs:
%   F (2*length(xdata)) function values at coeff
%   J (2*length(xdata) x length(coeff)) Jacobian matrix at coeff

    K = length(xdata); % number of angular frequencies
    Nxi = length(coeff)/2; % number of poles
        % Function definition for optimization
    coeff2xi = @(coeff)coeff(1:Nxi);
    coeff2mu = @(coeff)coeff((Nxi+1):end);
    H_optim = @(s,coeff)get1stOrderUnitaryFilter(s,coeff2xi(coeff),'type','high-pass')*coeff2mu(coeff);
    F = [real(H_optim(1i*xdata(:),coeff));imag(H_optim(1i*xdata(:),coeff))];
if nargout > 1   % Two output arguments
        % Matrix: (1i*omega_i+xi_j)^-1 (size K X Nxi)
    H = get1stOrderUnitaryFilter(1i*xdata(:),coeff2xi(coeff),'type','low-pass');
        % Matrix: -mu_j/(1i*omega_i+xi_j)^2 (size K x Nxi)
    G = -(H.*H).*kron(ones(length(H),1),transpose(coeff2mu(coeff))); % a bit costly in memory
        % Multiplication matrix [1i*omega_i]_i,j (size K x Nxi)
    Mult = kron(1i*xdata(:),ones(1,Nxi));
        % Jacobian matrix
    J = zeros(2*K,2*Nxi);
            % Real part w.r.t. xi (K*Nxi block)
    J(1:K,1:Nxi) = real(Mult.*G);
            % Real part w.r.t. mu (K*Nxi block)
    J(1:K,Nxi+(1:Nxi)) = real(Mult.*H);
            % Imag part w.r.t. xi (K*Nxi block)
    J(K+(1:K),1:Nxi) = imag(Mult.*G);
            % Imag part w.r.t. mu (K*Nxi block)
    J(K+(1:K),Nxi+(1:Nxi)) = imag(Mult.*H);
end
end

function [F,J] = CostFunRefl(coeff,xdata)
%fun Function to be used in call to lsqcurvefit.
% Cost function Jpsi.
% Outputs:
%   F (2*length(xdata)) function values at coeff
%   J (2*length(xdata) x length(coeff)) Jacobian matrix at coeff

    K = length(xdata); % number of angular frequencies
    Nxi = length(coeff)/2; % number of poles
        % Function definition for optimization
    coeff2xi = @(coeff)coeff(1:Nxi);
    coeff2mu = @(coeff)coeff((Nxi+1):end);
    H_optim = @(s,coeff)z2beta(get1stOrderUnitaryFilter(s,coeff2xi(coeff),'type','low-pass')*coeff2mu(coeff));
    F = [real(H_optim(1i*xdata(:),coeff));imag(H_optim(1i*xdata(:),coeff))];
if nargout > 1   % Two output arguments
        % Matrix: (1i*omega_i+xi_j)^-1 (size K X Nxi)
    H = get1stOrderUnitaryFilter(1i*xdata(:),coeff2xi(coeff),'type','low-pass');
        % Matrix: -mu_j/(1i*omega_i+xi_j)^2 (size K x Nxi)
    G = -(H.*H).*kron(ones(length(H),1),transpose(coeff2mu(coeff))); % a bit costly in memory
        % Matrix: (1+sum_n mu_n/(1i*omega_i+xi_n)_i,j (size K x Nxi)
    R = 1+get1stOrderUnitaryFilter(1i*xdata(:),coeff2xi(coeff),'type','low-pass')*coeff2mu(coeff);
    R = kron(R,ones(1,Nxi));
        % Jacobian matrix
    J = zeros(2*K,2*Nxi);
            % Real part w.r.t. xi (K*Nxi block)
    J(1:K,1:Nxi) = real(2*G./(R.^2));
            % Real part w.r.t. mu (K*Nxi block)
    J(1:K,Nxi+(1:Nxi)) = real(2*H./(R.^2));
            % Imag part w.r.t. xi (K*Nxi block)
    J(K+(1:K),1:Nxi) = imag(2*G./(R.^2));
            % Imag part w.r.t. mu (K*Nxi block)
    J(K+(1:K),Nxi+(1:Nxi)) = imag(2*H./(R.^2));
end
end

function [F,J] = CostFunReflNormalized(coeff,xdata,h_an)
%fun Function to be used in call to lsqcurvefit.
% Cost function Jpsi, normalized.
% Outputs:
%   F (2*length(xdata)) function values at coeff
%   J (2*length(xdata) x length(coeff)) Jacobian matrix at coeff
%   h_an (fun_handle) Standard diffusive function being approximated (Laplace)

    K = length(xdata); % number of angular frequencies
    Nxi = length(coeff)/2; % number of poles
        % Function definition for optimization
    coeff2xi = @(coeff)coeff(1:Nxi);
    coeff2mu = @(coeff)coeff((Nxi+1):end);
    H_optim = @(s,coeff)z2beta(get1stOrderUnitaryFilter(s,coeff2xi(coeff),'type','low-pass')*coeff2mu(coeff));
    tmp = H_optim(1i*xdata(:),coeff)./z2beta(h_an(1i*xdata(:)));
    F = [real(tmp);imag(tmp)];
if nargout > 1   % Two output arguments
        % Matrix: (1i*omega_i+xi_j)^-1 (size K X Nxi)
    H = get1stOrderUnitaryFilter(1i*xdata(:),coeff2xi(coeff),'type','low-pass');
        % Matrix: -mu_j/(1i*omega_i+xi_j)^2 (size K x Nxi)
    G = -(H.*H).*kron(ones(length(H),1),transpose(coeff2mu(coeff))); % a bit costly in memory
        % Matrix: (1+sum_n mu_n/(1i*omega_i+xi_n)_i,j (size K x Nxi)
    R = 1+get1stOrderUnitaryFilter(1i*xdata(:),coeff2xi(coeff),'type','low-pass')*coeff2mu(coeff);
    R = kron(R,ones(1,Nxi));
            % Mult matrix [1/z2beta(H(1i*omega_i))]_i,j (size K x Nxi)
    Mult = kron(1./z2beta(h_an(1i*xdata(:))),ones(1,Nxi));
        % Jacobian matrix
    J = zeros(2*K,2*Nxi);
            % Real part w.r.t. xi (K*Nxi block)
    J(1:K,1:Nxi) = real(Mult.*(2*G./(R.^2)));
            % Real part w.r.t. mu (K*Nxi block)
    J(1:K,Nxi+(1:Nxi)) = real(Mult.*(2*H./(R.^2)));
            % Imag part w.r.t. xi (K*Nxi block)
    J(K+(1:K),1:Nxi) = imag(Mult.*(2*G./(R.^2)));
            % Imag part w.r.t. mu (K*Nxi block)
    J(K+(1:K),Nxi+(1:Nxi)) = imag(Mult.*(2*H./(R.^2)));
end
end

function beta = z2beta(z)
    beta = (z-1)./(z+1);
end
