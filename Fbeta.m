%----------------------------------------------------------------------------------------
% File: Fbeta.m
%
% Goal: Compute Fbeta using the Nystrom interpolant
%
% Use: [fbeta] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin)
%
% Input: kind - see VieSolve.m file
%        T0 - step for iteration
%        t - row array of the evaluation points
%        a - solution of the linear system
%        x - row array of the zeros of the m-th Laguerre polynomial
%        w - row array of the corresponding Christoffel numbers
%        w1 - weights of the M-point Gauss-Laguerre rule with with M = 2048
%        E1 - array exp(-x1), with x1 array of the zeros of the M-th Laguerre polynomial
%        U - u(x), i.e weight function u computed ad the Laguerre zeros x
%        c - matrix (c_k(t_i))_{i = 1,...,length(t), k = 1,...,j}
%        type_OB - see VieSolve.m
%        varargin - see VieSolve.m
%
% Output: Fbeta - array of the solution of the VIE at the evaluation points t
%
% Recalls: IterativeNystromInterp.m
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
function [fbeta] = Fbeta(kind,T0,t,a,x,w,w1,E1,U,type_OB,varargin)

mm = length(t);
fbeta = zeros(mm,1);
% determine the number iter of iterations required for the computation of 
% fbeta at any evaluation time t(i), i=1:mm. 
% The iterative method is applied when t(i) > T0
for i = 1:mm
    if t(i) <= T0
        iter = 0;
    else
        v = (t(i)-mod(t(i),T0))/T0-1;
        if mod(t(i),T0) == 0
            iter = v-1;
        else
            iter = v;
        end
    end
    y = t(i)-iter*T0;
    % Compute fbeta using the Nystrom interpolant
    [fbeta(i)] = IterativeNystromInterp(kind,T0,iter,y,a,x,w,w1,E1,U,type_OB,varargin{:});
end