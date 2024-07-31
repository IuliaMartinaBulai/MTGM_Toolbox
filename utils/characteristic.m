%--------------------------------------------------------------------------
% File: characteristic.m
%
% Goal: define the characteristic function
%
% Use: [y] = characteristic(x,vbar)
%
% Input: x - row vector of values
%        vbar - lower bound of the volume of the metastases whose
%               cumulative number Nv is computed
%
% Output: [y] - characteristic function
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
function [y] = characteristic(x,vbar)
n = length(x);
y = zeros(n,1);
for i = 1:n
    if (x(i) < vbar)
        y(i) = 0;
    else
        y(i) = 1;
    end
end
