function H = get1stOrderUnitaryFilter(s,xi,varargin)
%get1stOrderUnitaryFilter 1/(s+xi) (low-pass) or s/(s+xi) (high-pass).
% Inputs:
%   s (N) Laplace variable (complex)
%   xi (P) pulsation (rad/s)
% Inputs (optional):
%   type (str) 'low-pass' (default) or 'high-pass'
% Output:
%   H (NxP)
%
% Copyright (c) 2018 Florian Monteghetti.
% This work is licensed under the terms of the MIT license.  
% For a copy, see <https://opensource.org/licenses/MIT>.


%   Slow version
%     for j=1:length(s)
%         H(j,:) = (s(j)+xi).^(-1); or H(j,:) = s(j)*(s(j)+xi).^(-1);
%     end

        % Optional argument
    p = inputParser; p.FunctionName = mfilename;
    addParameter(p,'type','low-pass',@(x)(validateattributes(x,{'char'},{})));
    parse(p,varargin{:});  oA = p.Results;
    type = oA.type;


    if isempty(s) || isempty(xi)
        H = [];
    else
        validateattributes(s,{'numeric'},{'vector'},1);
        validateattributes(xi,{'numeric'},{'vector'},2);
            % Convert to row vector (if necessary)
        s=s(:);
        xi=xi(:);

        H = kron(s,ones(1,length(xi))) + kron(ones(length(s),1),transpose(xi));
        if strcmp(type,'high-pass')
            H = kron(s,ones(1,length(xi)))./H;
        else
            H = 1./H;
        end
    end
end
