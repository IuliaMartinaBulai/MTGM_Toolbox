%----------------------------------------------------------------------------------------
% File: SystemBuild.m
%
% Goal: build the matrix and the right-hand side terms of the linear systems
%
% Use: [A,b] = SystemBuild(kind,T0,iter,a,x,w,w1,E1,U,c,type_OB,varargin)
%
% Input: kind - see VieSolve.m file
%        T0 - step for iteration
%        iter - number of iterations
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
% Output: A - matrix of the linear system
%         b - right-hand side of the linear system
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
function [A,b] = SystemBuild(kind,T0,iter,a,x,w,w1,E1,U,c,type_OB,varargin)
j = length(U);
xn = x+iter*T0;
Y = xn(1:j).*E1';
% build the matrix A of the linear system
KK = K(kind,xn(1:j),xn(1:j)',varargin{:});
A = eye(j)-KK.*c.*U';
% build the right-hand side term b of the linear system for Fbeta
H = 0;
cc = cK(T0,j,x,w,U)';
for i = 0:iter-1
    xi = x+i*T0;
    Kmat = K(kind,xi(1:j),xn(1:j)',varargin{:});
    H = H+Kmat*(cc.*a(i*j+1:(i+1)*j));
end

b = ((Gbeta(kind,xn(1:j),w1,Y,type_OB,varargin{:})+H').*U)';
end
