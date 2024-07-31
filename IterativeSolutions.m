%----------------------------------------------------------------------------------------
% File: IterativeSolutions.m
%
% Goal: Compute the solutions of the linear systems and their condition numbers 
%       at the iteration step iter+1
%
% Use: [a,C] = IterativeSolutions(kind,T0,T,m,x,w,w1,E1,U,type_OB,varargin)
%
% Input: kind -  see VieSolve.m file
%        T0 - step for iteration
%        T - greatest evaluation point
%        m - number of knots
%        x - row array of the zeros of the m-th Laguerre polynomial
%        w - row array of the corresponding Christoffel numbers
%        w1 - weights of the M-point Gauss-Laguerre rule with with M = 2048
%        E1 - array exp(-x1), with x1 array of the zeros of the M-th Laguerre polynomial
%        U - u(x), i.e weight function u computed ad the Laguerre zeros x
%        type_OB - see VieSolve.m
%        varargin - see VieSolve.m
%
% Output: a - column array of length 2*j*(iter+1) whose entries are the
%             solutions of the linear systems solved at the iteration step
%             iter+1
%         C - condition number of the last solved linear system
%
% Recalls: NystromMethod.m
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
function [a,C] = IterativeSolutions(kind,T0,T,m,x,w,w1,E1,U,type_OB,varargin)
if T <= T0
    [a,C] = NystromMethod(kind,0,0,[],m,x,w,w1,E1,U,type_OB,varargin{:});
else 
    disp('Number of iterations');
    iter = (T-mod(T,T0))/T0-1
    [a] = NystromMethod(kind,0,0,[],m,x,w,w1,E1,U,type_OB,varargin{:});
    for i = 1:iter  
        disp('Iteration step');
        i
        [a,C] = NystromMethod(kind,T0,i,a,m,x,w,w1,E1,U,type_OB,varargin{:});    
    end
    if mod(T,T0) ~= 0
        [a,C] = NystromMethod(kind,T0,(iter+1),a,m,x,w,w1,E1,U,type_OB,varargin{:});
    end
end
