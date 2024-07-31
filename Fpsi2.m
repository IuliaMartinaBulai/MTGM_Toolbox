%--------------------------------------------------------------------------
% File: Fpsi.m
%
% Goal: Compute the total metastatic mass, M, and the cumulative number of 
%       metastases, N, for the 2D non-autonomous case
%
% Use: [M,N] = Fpsi(kind,T0,t,a,x,w,w1,E,E1,U,type_OB,varargin)
%
% Input: kind - can assume values:
%             - 0 - to solve a Volterra Integral Equation
%             - 1 - to compute the cumulative number of metastases and
%                   the total metastatic mass for the 1D metastatic tumor
%                   growth model
%             - 2 - to compute the cumulative number of metastases and the
%                   total metastatic mass for the 1D metastatic tumor
%                   growth model solving the 1D ODE
%             - 3 - to compute the cumulative number of metastases and the
%                   total metastatic mass for the 2D metastatic tumor
%                   growth model with treatment solving the 2D ODE
%                   non-autonomous system OR
%             - 4 - to compute the cumulative number of metastases and the
%                   total metastatic mass for the 2D metastatic tumor
%                   growth model solving the 2D ODE autonomous system.
%        T0 - step time for iteration
%        t -  row array of the evaluation times expressed in days
%        a - solution of the linear system
%        x - Laguerre zeros
%        w - Christoffel numbers
%        w1 - weights of the N-point Gauss-Laguerre rule with N points
%        E - array exp(-x1) with x1 array of N Laguerre zeros
%        E1 - array exp(-x1) with x1 array of N1 Laguerre zeros
%        U - u(x) weight function u computed ad the Laguerre zeros x
%        type_OB - can assume value:
%                - 1 - to compute the volume of the metastatic mass OR
%                - 2 - to compute the cumulative number of metastases
%        varargin - has different elements depending on the value of kind,
%                   for kind = 0 see vie.m file
%                   for kind = 1 see tumorGrowthFunctions.m
%                   for kind = 2 see tumorGrowthFunctionsOdeSolve.m
%                   for kind = 3 and kind = 4 see tumorTreatmentParam.m
%
% Output: M -  1xT array of the approximate values of the metastatic mass
%              at t = [1:T] days
%         N - 1xT array of the approximate values of the cumulative 
%              number of metastases whose volume is larger than Vbar at
%              t = [1:T] days
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
function [M,N] = Fpsi2(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin)
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
        i=1;
        add1=1;
        add2=1;
        while (i< length(w)) && (abs(add1)>0.5e-25) && (abs(add2)>0.5e-25)
            [fbeta] = Fbeta(kind,T0,y*E1(i),a,x,w,w1,E1,U,type_OB,varargin{:})
            [~,sol] = ode23t(@(tt,yy)odefun2Dm (tt, yy, param),[y*E1(i) y],y0_m,options);
            [~,sol2] = ode23t(@(tt,yy)odefun2Dp (tt, yy, param),[0  y*E1(i)],y0_p,options);
            add1=w1(i)*sol(end,1)*(be_p(sol2(end,1))+fbeta)
            M(ss) = M(ss)+add1;
            add2= w1(i)*characteristic(sol(end,1),vbar)*(be_p(sol2(end,1))+fbeta)
            N(ss) = N(ss)+add2;
        end
    end
    M(ss) = y*M(ss);
    N(ss) = y*N(ss);
end
M = M';
N = N';
end