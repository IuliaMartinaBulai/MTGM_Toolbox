%----------------------------------------------------------------------------------------
% File: tumorTreatmentParam.m
%
% Goal: Create a database of possible tumor treatments
%
% Use: [g,y0_p,y0_m,param,vbar,options] = tumorTreatmentParam(varargin)
%
% Input: varargin - 1X12 cell containing 6 couples:
%                   'treatment_type', value1 - no_treat or type of treatment
%                   'emission_p',mu_p - colonization coefficient of the primary tumor
%                   'emission_m',mu_m  - colonization coefficient of the metastases
%                   'y0_p', [x0_p,theta0_p] - initial condition of the size of the tumor
%                                             (x0_p) and of the carrying capacity
%                                             (theta0_p) of the primary tumor
%                   'y0_m', [x0_m,theta0_m] - initial condition of the size of the tumor
%                                             (x0_m) and of the carrying capacity
%                                             (theta0_m) of the secondary tumor
%                   'vbar', vbar - lower bound of the volume of the metastases whose
%                                  cumulative number N is computed
%
% Output: g - 1X4 cell with
%             g{1} emission rate function for the primary tumor be_p
%             g{2} emission rate function for the metastases be_m
%             g{3} weight function u of the weighted space
%         y0_p - [x0_p,theta0_p] initial conditions of the size of the tumor (x0_p) and
%                of the variable carrying capacity (theta0_p) for the primary tumor
%         y0_m - [x0_m,theta0_m] initial conditions of the size of the tumor (x0_m) and
%                 of the variable carrying capacity (theta0_m) for the secondary tumor
%         param - [DA, clrA, DC, clrC, a, cc, d, e, h, xmin, thetamin, timesA, timesC]
%                  parameters related to the treatment type
%         vbar - lower bound of the volume of the metastases whose cumulative number N
%                is computed
%
% Recalls: argselectAssign.m, argselectCheck.m
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
function [g,y0_p,y0_m,param,vbar,options] = tumorTreatmentParam(varargin)

% default parameter values
alpha = 2/3;
mu_p = 10^-3; % day^(-1)mm^(-3 al)
mu_m = 10^-3; % day^(-1)mm^(-3 al)
% default case study
value1 = 'no_treat';
% default initial values
x0_p = 200; 
x0_m = 10^-6;
theta0_p = 625;
y0_p = [x0_p, theta0_p];
y0_m = [x0_m, theta0_p];
vbar = 0;
options = [];
%options=odeset('Reltol',1e-10,'Abstol',1e-10);

control_params = {'treatment_type', value1,'emission_p',mu_p,'emission_m',mu_m,...
    'y0_p', [x0_p, theta0_p],'y0_m', [x0_m, theta0_p], 'vbar', vbar};
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

% default parameter values (releated to the treatment)
a = 0.192;
cc = 5.85;
d = 8.73*10^-3;
e = 0;
h = 0;
timesA = (5:1:10);
timesC = (10:1:15);
DC = 1;
thetamin = 10^-6;
xmin = 10^-6;
DA = 20;
clrC = 1;
clrA = 1.7;

switch treatment_type
    case 'no_treat'
    case 'end_1'
        e = 0.66;
    case 'end_2'
        e = 0.66;
        timesA = (5:2:15);
    case 'end_3'
        e = 0.66;
        timesA=(5:0.5:7.5);
    case 'end_4'
        e = 0.66;
        DA = 5;
    case 'end_5'
        e = 0.66;
        DA = 10;
    case 'end_6'
        e = 0.66;
        DA = 30;
    case 'end_7'
        e = 0.66;
        DA = 40;
    case 'angiostatin'
        e = 0.15;
        clrA = 0.38;
    case 'TNP-470'
        e = 1.3;
        clrA = 10.1;
        DA = 30;
        timesA = (4:2:10);
    case 'AACT'
        e = 0.66;
        h = 1;
    case 'CTAA'
        e = 0.66;
        timesA = (10:1:15);
        h = 1;
        timesC = (5:1:10);
end
param = [DA, clrA, DC, clrC, a, cc, d, e, h, xmin, thetamin, timesA, timesC];
end
