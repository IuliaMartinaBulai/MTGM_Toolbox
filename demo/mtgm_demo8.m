%--------------------------------------------------------------------------
% mtgm_demo8:
%
% This demo solves the VIE in order to compare the effect of the combination 
% of two different drugs, an anti-angiogenic (AA) one and a cytotoxic one 
% (CT) on the cumulative number of metastases, N(t), and on the total 
% metastatic mass, M(t).
%
% This reproduces figure 5 from "Numerical solution of metastatic tumor 
% growth models with treatment", IM Bulai, MC De Bonis, C Laurita, 2024.
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

fprintf('Welcome to MTGM demo #8\n');
fprintf(['Compute the metastatic mass, M(t), and the cumulative number of ' ...
    'metastases, N(t) \n'])
fprintf('2D non-autonomous model with treatment \n')
addpath('..');
addpath('../utils');

% time in days
T = 15;
T0 = 10;
t = (1:T);
k = 6;
kind = 3;
type_OB = 0;

% No treatment 
tic
fprintf('No treatment \n');
varargin_no_treat = {};
[M_no_treat,N_no_treat,j_no_treat,C_no_treat] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_no_treat{:});
toc
fprintf(' Index j and condition number C with No treatment \n');
[j_no_treat, C_no_treat]

% AACT
fprintf('AACT \n');
varargin_end_AACT = {'treatment_type','AACT'};
tic
[M_end_AACT,N_end_AACT,j_end_AACT,C_end_AACT] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_AACT{:});
toc
fprintf(' Index j and condition number C with AACT \n');
[j_end_AACT, C_end_AACT]

% CTAA
fprintf('CTAA \n');
varargin_CTAA = {'treatment_type','CTAA'};
tic
[M_CTAA,N_CTAA,j_CTAA,C_CTAA] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_CTAA{:});
toc
fprintf(' Index j and condition number C with CTAA \n');
[j_CTAA, C_CTAA]

% change the directory to save the figures: 
cd Figures

% plot and save the results
fprintf('Plot metastatic mass \n');
figure
plot(t,M_no_treat,'b-*', t,M_CTAA,'r-+',t,M_end_AACT,'k-o','linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
legend('Without treatment','CT + AA','AA + CT','Location','best');
saveas(gcf, 'Met_mass_demo8', 'fig');
saveas(gcf, 'Met_mass_demo8', 'epsc');

fprintf('Plot metastatic mass zoom \n');
figure
plot(t,M_CTAA,'r-+',t,M_end_AACT,'k-o','linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
legend('CT + AA','AA + CT','Location','best');
saveas(gcf, 'Met_mass_demo8_zoom', 'fig');
saveas(gcf, 'Met_mass_demo8_zoom', 'epsc');

fprintf('Plot cumulative number of metastases \n');
figure
plot(t,N_no_treat,'b-*',t, N_CTAA,'r-+',t,N_end_AACT,'k-o','linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Cumulative number')
legend('Without treatment','CT + AA','AA + CT','Location','best');
saveas(gcf, 'Cum_numb_demo8', 'fig');
saveas(gcf, 'Cum_numb_demo8', 'epsc');