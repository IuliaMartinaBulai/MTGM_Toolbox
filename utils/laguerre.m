%--------------------------------------------------------------------------
% File: laguerre.m
%
% Goal: compute the sequence of the orthonormal Laguerre polynomials using the
%       recurrence formula 
%             b p (x) = (x - a ) p   (x) - b   p   (x),    j=1,2,...,n,
%              j j            j   j-1       j-1 j-2
%
% Use: [lag2] = laguerre(n,alpha,x)
%
% Input: n - degree of the last polynomial of the sequence
%                                                                         alpha  -x                                                                                                                                                      
%        alpha - parameter > -1 defining the Laguerre weight function w  (x)= x      e   
%                                                                     alpha
%        x - evaluation point
%
% Output: lag2 - orthonormal Laguerre polynomials
%
% Recalls: class1.m
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
function [lag2] = laguerre(n,alpha,x)
[a,b,muzero] = class1(6, n+1, alpha, 0);
lag0 = (muzero)^(-0.5); 
lag(1) = ((x-a(1))*lag0)/b(1);
lag(2) = ((x-a(2))*lag(1)-b(1)*lag0)/b(2);
for j = 3:n
   lag(j) = ((x-a(j))*lag(j-1)-b(j-1)*lag(j-2))/b(j);
end 
lag2 = [lag0,lag];
end