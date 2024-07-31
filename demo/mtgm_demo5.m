%--------------------------------------------------------------------------
% mtgm_demo5:
%
% This demo solves the VIE in order to compute the cumulative number of
% metastases, N(t), and the total metastatic mass M(t) for a metastatic 
% tumor growth model. Differently than for demo3 here the 2D autonomous
% model is used.
%
% This reproduces figure ... from "A new MATLAB software for numerical
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
fprintf('Welcome to MTGM demo #5\n');
fprintf(['Compute the metastatic mass, M(t), and the cumulative number of ' ...
    'metastases, N(t) \n'])
fprintf('2D autonomous model without treatment \n')
addpath('..');
addpath('../utils');

T = 40;
T0 = 20;
t = (1:T);
k = 6;
kind = 4;
varargin = {'y0_p', [10^-6, 625]};

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
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
saveas(gcf, 'Met_mass_demo5', 'fig');
saveas(gcf, 'Met_mass_demo5', 'epsc');


fprintf('Plot cumulative number of metastases \n');
figure
plot(t,N, 'linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Cumulative number')
saveas(gcf, 'Cum_numb_demo5', 'fig');
saveas(gcf, 'Cum_numb_demo5', 'epsc');

