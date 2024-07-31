%----------------------------------------------------------------------------------------
% File: NystromMethod.m
%
% Goal: Compute the solution of the linear system and its condition number 
%       at the iteration step iter+1 and update the array a of the
%       solutions
%
% Use: [a,C] = NystromMethod(kind,T0,iter,a,m,x,w,w1,E1,U,type_OB,varargin)
%
% Input: kind - see VieSolve.m file
%        T0 - step for iteration
%        iter - number of iterations
%        a - column array of length 2*j*iter whose entries are the
%            solutions of the linear systems solved at the iteration step
%            iter
%        m - number of knots
%        x - row array of the zeros of the m-th Laguerre polynomial
%        w - row array of the corresponding Christoffel numbers
%        w1 - weights of the M-point Gauss-Laguerre rule with with M = 2048
%        E1 - array exp(-x1), with x1 array of the zeros of the M-th Laguerre polynomial
%        U - u(x), i.e weight function u computed ad the Laguerre zeros x
%        type_OB - see VieSolve.m
%        varargin - see VieSolve.m
%
% Output: a - column array of length 2*j*(iter+1) whose entries are all the
%             solutions of the linear systems solved at the iteration step
%             iter+1
%         C - condition number of the linear system at the iteration step
%             iter+1
%
% Recalls: cK.m, SystemBuild.m
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

function [a,C] = NystromMethod(kind,T0,iter,a,m,x,w,w1,E1,U,type_OB,varargin)
theta = 0.25;
j = length(U);
% for n<=1024 upload c (previously computed and saved as .mat) otherwise
% compute it using the function cK.m
if (m==64 && theta==0.25)
    c = importdata('c64.mat');
elseif (m==128 && theta==0.25)
    c = importdata('c128.mat');
elseif (m==256 && theta==0.25)
    c = importdata('c256.mat');
elseif (m==512 && theta==0.25)
    c = importdata('c512.mat');
elseif (m==1024 && theta==0.25)
    c = importdata('c1024.mat');
else
    [c] = cK(x(1:j),j,x,w,U);
end

% build the matrix and the right-hand side terms of the linear systems for
% the total metastatic mass and the cumulative number of metastases
[A,b] = SystemBuild(kind,T0,iter,a,x,w,w1,E1,U,c,type_OB,varargin{:});
% compute the condition number of the matrix A
C = cond(A,inf);
% compute the solution of the linear systems
asol = A\b;
a = [a; asol];
end