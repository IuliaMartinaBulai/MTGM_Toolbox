%----------------------------------------------------------------------------------------
% File: cK.m
%
% Goal: Exact computation of the integrals
%                t
%    c_k(t) =   | l_{n+1,k}(z) dz,  k = 1,...,j
%               0
%
% Use: [c] = cK(t,j,x,w,U)
%
% Input: t - row array of evaluation points
%        j - size of the solved linear system
%        x - row array of the zeros of the m-th Laguerre polynomial
%        w - row array of the corresponding Christoffel numbers
%        U - u(x), i.e weight function u computed at the Laguerre zeros x
%
% Output: c - matrix (c_k(t_i))_{i = 1,...,length(t), k = 1,...,j}
%
% Recalls: class1.m, laguerre.m
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
function [c] = cK(t,j,x,w,U)

n = length(x);
[a,b,~] = class1(6,n+1,0,0);
m = length(t);
c = zeros(m,j);
mom = zeros(n,m);
momento = zeros(1,n-1);

for nu = 1:m
    lag = laguerre(n+1,0,t(nu));
    mom0 = (4*n-a(1))*lag(1)+(4*n-a(1)-b(1))*lag(2)-b(1)*lag(3);
    for i = 1:n-1
        momento(i) = (4*n-a(i+1)-b(i))*lag(i+1)+(4*n-a(i+1)-b(i+1))*lag(i+2)...
            -b(i)*lag(i)-b(i+1)*lag(i+3);
    end
    mom(:,nu) = [mom0,momento]';

    for k = 1:j
        lagzk = laguerre(n-1,0,x(k));
        c(nu,k) = w(k)*lagzk*mom(:,nu)/(U(k)*(4*n-x(k)));
    end
end
