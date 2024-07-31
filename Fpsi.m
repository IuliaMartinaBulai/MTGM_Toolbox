%----------------------------------------------------------------------------------------
% File: Fpsi.m
%
% Goal: Compute the total metastatic mass, M, and the cumulative number of metastases,
%       N, for the 2D non-autonomous case
%
% Use: [M,N] = Fpsi(T0,t,a,x,w,w1,E,E1,U,varargin)
%
% Input: T0 - step for iteration
%        t - row array of the evaluation points
%        a - solution of the linear system
%        x - row array of the zeros of the m-th Laguerre polynomial
%        w - row array of the corresponding Christoffel numbers
%        w1 - weights of the M-point Gauss-Laguerre rule with with M = 2048
%        E - array exp(-x)
%        E1 - array exp(-x1), with x1 array of the zeros of the M-th Laguerre polynomial
%        U - u(x), i.e weight function u computed ad the Laguerre zeros x
%        varargin - see VieSolve.m
%
% Output: M -  1xT array of the approximate values of the metastatic mass at t = [1:T]
%              days
%         N - 1xT array of the approximate values of the cumulative number of metastases
%              whose volume is larger than Vbar at  t = [1:T] days
%
% Recalls: tumorTreatmentParam.m, Fbeta.m, ode23t.m, odefun2Dm .m,
%          odefun2Dp.m, characteristic.m
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
function [M,N] = Fpsi(T0,t,a,x,w,w1,E,E1,U,varargin)
kind = 3;
type_OB = 0;
j = length(U);
[g,y0_p,y0_m,param,vbar,options] = tumorTreatmentParam(varargin{:});
be_p = g{1};
mm = length(t);
% compute the volume of the metastatic mass and the cumulative number of
% metastases greater than vbar
M = zeros(1,mm);
N = zeros(1,mm);
parfor ss = 1:mm
    y = t(ss)
    if y ~= 0
        for i = 1:j
            [fbeta] = Fbeta(kind,T0,y*E(i),a,x,w,w1,E1,U,type_OB,varargin{:})
            [~,sol] = ode23t(@(tt,yy)odefun2Dm (tt, yy, param),[y*E(i) y],y0_m,options);
            [~,sol2] = ode23t(@(tt,yy)odefun2Dp (tt, yy, param),[0  y*E(i)],y0_p,options);
            M(ss) = M(ss)+w(i)*sol(end,1)*(be_p(sol2(end,1))+fbeta);
            N(ss) = N(ss)+w(i)*characteristic(sol(end,1),vbar)*(be_p(sol2(end,1))+fbeta);
        end
    end
    M(ss) = y*M(ss);
    N(ss) = y*N(ss);
end
M = M';
N = N';
end