%--------------------------------------------------------------------------
% mtgm_demo2:  
%
% This demo solves the VIE in order to compute the metastatic mass, M(t), 
% and the cumulative number of metastases, N(t), for breast tumor data
% assuming that both the primary and secondary tumors growth laws are 
% Gompertz laws. Differently from mtgm_demo1 here the ODE equation is 
% solved numerically. 
%
% This reproduces figure 4 from "A new MATLAB software for numerical
% computation of biological observables for metastatic tumor growth",
% IM Bulai, MC De Bonis, C Laurita, 2024.
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
clc
clear 
close all
fprintf('Welcome to MTGM demo #2\n');
addpath('..');
addpath('../utils');

fprintf(['compute the metastatic mass, M(t), and the cumulative number of ' ...
    'metastases, N(t) \n']);
fprintf('Case study: breast \n');
T = 60;
T0 = 30;
t = (1:T);
k = 8;
kind = 2;
varargin = {};

tic
type_OB = 1;
fprintf('Compute the metastatic mass M \n');
[M,~,j_M,C_M] = VieSolve(kind,T0,t,2^k,type_OB,varargin{:});
type_OB = 2;
fprintf('Compute the cumulative number of metastases N \n');
[~,N,j_N,C_N] = VieSolve(kind,T0,t,2^k,type_OB,varargin{:});
toc

fprintf('Index j and condition number C for the metastatic mass M \n');
[j_M, C_M]

fprintf('Index j and condition number C for the cumulative number of metastases \n');
[j_N, C_N]

% change the directory to save the figures: 
cd Figures

% plot and save the results
fprintf('Plot metastatic mass \n');
figure
plot(t,M, 'linewidth',2)
set(gca,'fontsize',16)
title('Breast')
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
saveas(gcf, 'Met_mass_demo2', 'fig');
saveas(gcf, 'Met_mass_demo2', 'epsc');

fprintf('Plot cumulative number of metastases \n');
figure
plot(t,N, 'linewidth',2)
set(gca,'fontsize',16)
title('Breast')
xlabel('Time (days)')
ylabel('Cumulative number')
saveas(gcf, 'Cum_numb_demo2', 'fig');
saveas(gcf, 'Cum_numb_demo2', 'epsc');