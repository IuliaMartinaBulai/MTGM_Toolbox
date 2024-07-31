%----------------------------------------------------------------------------------------
% File: IterativeNystromInterp.m
%
% Goal: Compute the Nystrom interpolants, solutions of VIEs, for Fbeta at the iteration step iter
%
% Use: [Fbeta] = IterativeNystromInterp(kind,T0,iter,t,a,x,w,w1,E,U,...
%                type_OB,varargin)
%
% Input: kind - see VieSolve.m file
%        T0 - step for iteration
%        iter - number of iterations
%        t - row array of the evaluation points
%        a - solution of the linear system
%        x - row array of the zeros of the m-th Laguerre polynomial
%        w - row array of the corresponding Christoffel numbers
%        w1 - weights of the M-point Gauss-Laguerre rule with with M = 2048
%        E - array exp(-x)
%        U - u(x), i.e weight function u computed ad the Laguerre zeros x
%        type_OB - see VieSolve.m
%        varargin - see VieSolve.m
%
% Output: Fbeta - array of the solution of the VIE at the evaluation points t at
%                 the iteration step iter
%
% Recalls: K.m, cK.m, Gbeta.m
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
function [Fbeta] = IterativeNystromInterp(kind,T0,iter,t,a,x,w,w1,E,U,type_OB,varargin)
j = length(U);
xn = x+iter*T0;
tn = t+iter*T0;
Y = tn.*E';
% compute the modified moments
KK = K(kind,xn(1:j),tn',varargin{:});
C = KK.*cK(t,j,x,w,U);
H = 0;
cc = cK(T0,j,x,w,U)';
for i = 0:iter-1
    xi = x+i*T0;
    Kmat = K(kind,xi(1:j),tn',varargin{:});
    H = H+Kmat*(cc.*a(i*j+1:(i+1)*j));
end
asol = a(iter*j+1:(iter+1)*j);
% compute the Nystrom interpolant for Fbeta
Fbeta = Gbeta(kind,tn,w1,Y,type_OB,varargin{:})'+H+C*asol;
end