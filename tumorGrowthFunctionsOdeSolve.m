%----------------------------------------------------------------------------------------
% File: tumorGrowthFunctionsOdeSolve.m
%
% Goal: Data for the primary and metastatic tumor growth
%
% Use: [g,y0_m,y0_p,param,vbar,options] = tumorGrowthFunctionsOdeSolve(varargin)
%
% Input: varargin - 1X8 cell containing 4 couples:
%                   'emission_p',mu_p - colonization coefficient of the primary tumor
%                   'emission_m',mu_m - colonization coefficient of the metastases
%                   'v_p0', vp0 - volume of the primary tumor at time t = 0
%                   'v_m0', vm0 - volume of the newly created metastases
%
% Output: g - 1X3 cell with:
%             g{1} emission rate function for the primary tumor be_p
%             g{2} emission rate function for the metastases be_m
%             g{3} weight function u of the weighted space
%         y0_p - [x0_p,theta0_p] initial conditions of the size of the 
%                tumor (x0_p) and of the variable carrying capacity 
%                (theta0_p) for the primary tumor
%         y0_m - [x0_m,theta0_m] initial conditions of the size of the 
%                tumor (x0_m) and of the variable carrying capacity 
%                (theta0_m) for the secondary tumor
%         param - [DA, clrA, DC, clrC, a, cc, d, e, h, xmin, thetamin, 
%                 timesA, timesC] parameters related to the treatment type
%         vbar - lower bound of the volume of the metastases whose 
%                cumulative number N is computed
%         options - integration settings for the ODE solver
%
% Recalls: argselectAssign.m, argselectCheck.m
%
% Authors: IM Bulai, MC De Bonis, C Laurita
% Date last modified: January, 2025
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
function [g,y0_m,y0_p,param,vbar,options] = tumorGrowthFunctionsOdeSolve(varargin)

% default parameter values
alpha = 2/3;
mu_p = 10^-3; % day^(-1)mm^(-3 al)
mu_m = 10^-3; % day^(-1)mm^(-3 al)

% default initial values
y0_m = 10^(-6); % mm^3
y0_p = 1;
vbar = 10^(-6);

% parameters for the Gompertz law
a = 0.56;
beta = 0.0719;

options = [];

control_params = {'emission_p',mu_p,'emission_m',mu_m,...
    'y0_m', y0_m,'y0_p', y0_p, 'vbar', vbar};
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

% emission rate function for the primary tumor be_p
g{1} = @(t) emission_p*t.^alpha;
% emission rate function for the metastases be_m
g{2} = @(t) emission_m*t.^alpha;
% u(s) in the paper, the weight function of the weighted space where
% the solutions of the VIEs lives
g{3} = @(y) y.^(1/4).*exp(-y/2);

param = [a, beta];
end
