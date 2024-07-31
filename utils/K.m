%----------------------------------------------------------------------------------------
% File: K.m
%
% Goal: compute the kernel of the VIE
%
% Use: [Kt] = K(kind,x1,x2,varargin)
%
% Input: kind - see VieSolve.m file
%        x1 - row array of the evaluation points
%        x2 - row array of the evaluation points
%        varargin - see VieSolve.m
%
% Output: Kt - matrix K(x1(1:j),x2(1:j)')
%
% Recalls: vie.m, tumorGrowthFunctions.m, tumorGrowthFunctionsOdeSolve.m, ode23t.m,
%          odefun1Dm.m, tumorTreatmentParam.m, odefun2Dm.m
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

function [Kt] = K(kind,x1,x2,varargin)

% kind = 0: Volterra Integral Equation
if kind == 0
    [g] = vie();
    s = x1;
    t = x2;
g1 = g{1};
Kt = g1(s,t);
    % kind = 1: Cumulative number of metastases and total metastatic mass for
    %           the 1D metastatic tumor growth model
elseif kind == 1
    [g,gamma,be_m,vm] = tumorGrowthFunctions(varargin{:});
    % k(s,t) in the paper, the kernel of the VIE
    s = x1;
    t = x2;
    Kt = be_m(vm(t-s));
    % kind = 2: Cumulative number of metastases and total metastatic mass
    %           for the 1D metastatic tumor growth model solving the 1D ODE
elseif kind == 2
    [g,y0_m,y0_p,param,options] = tumorGrowthFunctionsOdeSolve(varargin{:});
    be_p = g{1};
    be_m = g{2};
    m = length(x1);
    n = length(x2);
    Kt = zeros(n,m);
    parfor i = 1:n
        for k = 1:m
            if (x2(i)-x1(k)) >0
                [~,sol1]=ode23t(@(tt,yy)odefun1Dm (tt, yy, param),[0 x2(i)-x1(k)],y0_m,options);
                Kt(i,k) = be_m(sol1(end,1));
            elseif (x2(i)-x1(k)) ==0
                Kt(i,k) = be_m(y0_m);
            end
        end
    end
    % kind = 3: Cumulative number of metastases and total metastatic mass for
    %            the 2D metastatic tumor growth model with treatment solving
    %            the 2D ODE non-autonomous system
elseif kind == 3 
    [g,y0_p,y0_m,param,vbar,options] = tumorTreatmentParam(varargin{:});
    be_p = g{1};
    be_m = g{2};
    m = length(x1);
    n = length(x2);
    Kt = zeros(n,m);
    parfor i = 1:n
        for k = 1:m
            if x1(k) ~= x2(i) && x2(i)>x1(k)
                [~,sol1] = ode23t(@(tt,yy)odefun2Dm (tt,yy,param),[x1(k) x2(i)],y0_m,options);
                Kt(i,k) = be_m(sol1(end,1));
            end
        end
    end
    % kind = 4: Cumulative number of metastases and total metastatic mass for
    %           the 2D metastatic tumor growth model solving the 2D ODE
    %           autonomous system.
    elseif kind == 4
    [g,y0_p,y0_m,param,vbar,options] = tumorTreatmentParam(varargin{:});        
    be_p = g{1};
    be_m = g{2};
    m = length(x1);
    n = length(x2);
    Kt = zeros(n,m);
    parfor i = 1:n
        for k = 1:m
            if x1(k) ~= x2(i) && x2(i)>x1(k)
                [~,sol1]=ode23t(@(tt,yy)odefun2Dm (tt,yy,param),[x1(k) x2(i)],y0_m,options);
                Kt(i,k) = be_m(sol1(end,1));
            end
        end
    end
end
end





