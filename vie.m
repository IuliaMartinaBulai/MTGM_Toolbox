%----------------------------------------------------------------------------------------
% File: vie.m
%t
% Goal:  Define the analytical expressions of the kernel, the right-hand side term and
%        the weight function u of the weighted space, for the given VIE
%
% Use: [g] = vie()
%
% Input: -
% Output: g - 1X4 cell with
%             g{1} kernel K(s,t) of the VIE
%             g{2} right-hand side of the VIE
%             g{3} weight function u of the weighted space
%
% Recalls: -
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
function [g] = vie()
% kernel, right-hand side term and Laguerre weight
g{1} = @(s,t) exp (-s/2)./(s+t+10).^2; % K
g{2} = @(t) t-1/2*(-2+(exp(-t/2).*(10+t))./(5+t)+exp(5+t/2).*(12+t).*(real(-expint(5+t))-real(-expint((5+t/2))))); % G
g{3} = @(y) y.^(1/4).*exp(-y/2); % u
% exact solution f = @(t) t;
end