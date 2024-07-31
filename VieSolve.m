%----------------------------------------------------------------------------------------
% File: VieSolve.m
%
% Goal: Compute the required solution of the Volterra Integral Equation, the dimension
%       of the solved linear system and its condition number
%
% Use: [M,N,j,C] = VieSolve(kind,T0,t,m,type_OB,varargin)
%
% Input: kind - can assume values:
%             - 0 - to solve a Volterra Integral Equation
%             - 1 - to compute the cumulative number of metastases and the total metastatic
%                   mass for the 1D metastatic tumor growth model assuming that the growth
%                   laws for the primary and secondary tumor are known analytically;
%             - 2 - to compute the cumulative number of metastases and the total metastatic
%                   mass for the 1D metastatic tumor growth model solving numerically the
%                   differential problem describing the growth of the tumor;
%             - 3 - to compute the cumulative number of metastases and the total metastatic
%                   mass for the 2D metastatic tumor growth model with treatment solving
%                   the 2D ODE non-autonomous system;
%             - 4 - to compute the cumulative number of metastases and the total metastatic
%                   mass for the 2D metastatic tumor growth model solving the 2D ODE
%                   autonomous system
%        T0 - step for iteration
%        t - row array of the evaluation points
%        m - number of knots
%        type_OB - can assume value:
%                - 0 - for kind = 0 and kind = 3
%                - 1 - to compute the volume of the metastatic mass
%                - 2 - to compute the cumulative number of metastases
%        varargin - has different elements depending on the value of kind,
%                   for kind = 0 see vie.m file
%                   for kind = 1 see tumorGrowthFunctions.m
%                   for kind = 2 see tumorGrowthFunctionsOdeSolve.m
%                   for kind = 3 and kind = 4 see tumorTreatmentParam.m
%
% Output: M - 1xT array of the approximate values of:
%                  -the solution of the VIE at the evaluation array t for
%                   kind = 0;
%                  -the metastatic mass at t = [1:T] days for kind = 1,2,3,4
%         N - empty array for kind = 0;
%             1xT array of the approximate values of the cumulative number
%             of metastases whose volume is larger than Vbar at t = [1:T]
%             days for kind = 1,2,3,4
%         j - dimension of the solved linear system
%         C - condition number of the solved linear system
%
% Recalls: gaussq.m, IterativeSolutions.m, Fbeta.m, Fpsi.m
%
% Authors: IM Bulai, MC De Bonis, C Laurita
% Date last modified: July, 2024
%
% This file is part of the MTGM toolbox
% Copyright (C) 2024, IM Bulai, MC De Bonis, C Laurita.
%
% The MTGM toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by the
% Free Software Foundation.
%
% The MTGM toolbox is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the MTGM toolbox. If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

function [M,N,j,C] = VieSolve(kind,T0,t,m,type_OB,varargin)
theta = 0.25;

[x,w] = gaussq(6,m,0,0,0,[0,0]);
E = exp(-x);

j = 1;
while (j <= m) && (x(j) < 4*m*theta)
    j = j+1;
end
j = j-1;

% u(s) in the paper, the weight function of the weighted space where
% the solutions of the VIEs lives
u = @(y) y.^(1/4).*exp(-y/2);
U = u(x(1:j));

% compute N = 2048 Laguerre knots and weights
[x1,w1] = gaussq(6,2048,0,0,0,[0,0]);
E1 = exp(-x1);
if kind == 0 || 1 || 2 || 4
    T = max(t);
elseif kind == 3
    T = max(t)*E(1);
    %T = max(t)*E1(1);
end

[a,C] = IterativeSolutions(kind,T0,T,m,x,w,w1,E1,U,type_OB,varargin{:});

if kind == 0
    if type_OB == 0
        % compute fm u
        [M] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin{:}).*u(t)';
        N = [];
    end

    % compute M and N, only for the 2D non-autonomous case Fpsi must be called
    % otherwise use Fbeta
elseif kind == 1
    if type_OB == 1
        [M] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin{:});
        N = [];
    elseif type_OB == 2
        [N] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin{:});
        M = [];
    end
elseif kind == 2
    if type_OB == 1
        [M] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin{:});
        N = [];
    elseif type_OB == 2
        [N] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin{:});
        M = [];
    end
elseif kind == 3
    [M,N] = Fpsi(T0,t,a,x,w,w1,E,E1,U,varargin{:});
    %[M,N] = Fpsi2(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin{:})
elseif kind == 4
    if type_OB == 1
        [M] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin{:});
        N = [];
    elseif type_OB == 2
        [N] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin{:});
        M = [];
    end
end
end



