%--------------------------------------------------------------------------
% File: odefun2Dp.m
%
% Goal: function for the numerical integration of the ODE system for the
%       primary tumor
%
% Use: [dy] = odefun2Dp(t,y,par)
%
% Input: t -  row array of the evaluation times expressed in days
%        y - unknown of the ODE system
%        par - [DA, clrA, DC, clrC, a, cc, d, e, h, xmin, thetamin,
%                timesA, timesC] parameters related to the treatment type
%
% Output: [dy] - numerical solution of the ODE system
%
% Recalls: heaviside.m, subplus.m
%
% Authors: IM Bulai, MC De Bonis, C Laurita
% Date last modified: July, 2024
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

function [dy] = odefun2Dp(t,y,par)

if length(par)==21
    DA = par(1);
    clrA = par(2);
    DC = par(3);
    clrC = par(4);
    a = par(5);
    c = par(6);
    d = par(7);
    e = par(8);
    h = par(9);
    xmin = par(10);
    thetamin = par(11);
    timesA = par(12:15);
    timesC = par(16:21);
else
    DA = par(1);
    clrA = par(2);
    DC = par(3);
    clrC = par(4);
    a = par(5);
    c = par(6);
    d = par(7);
    e = par(8);
    h = par(9);
    xmin = par(10);
    thetamin = par(11);
    timesA = par(12:17);
    timesC = par(18:23);
end

dy = zeros(2,1);
gammaA = @(t) DA.*sum(exp(-clrA*(t-timesA)).*heaviside(t-timesA));
gammaC = @(t) DC.*sum(exp(-clrC*(t-timesC)).*heaviside(t-timesC));
dy(1) = a*y(1)*log(y(2)/y(1))-h*gammaC(t)*subplus(y(1)-xmin);
dy(2) = c*y(1)-d*y(2)*y(1)^(2/3)-e*gammaA(t)*subplus(y(2)-thetamin);
end