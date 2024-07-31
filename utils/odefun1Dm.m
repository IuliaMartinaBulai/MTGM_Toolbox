%--------------------------------------------------------------------------
% File: odefun1Dm.m
%
% Goal: function for the numerical integration of the ODE equation for the
%       secondary tumor, here we have used Gompertz growth law only for
%       sake of comparison
%
% Use: [dy] = odefun1Dm(t,y,par)
%
% Input: t - row array of the evaluation times expressed in days
%        y - unknown of the ODE system
%        par - [a, beta] parameters related to the model
%
% Output: [dy] - numerical solution of the ODE equation
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
function [dy] = odefun1Dm(t,y,par)
a = par(1);
beta = par (2);
dy = zeros(1,1);
dy(1) = a*exp(-beta*t)*y(1);
end