%--------------------------------------------------------------------------
% mtgm_demo6:
%
% This demo solves the VIE in order to compute the cumulative number of
% metastases, N(t), and the total metastatic mass M(t) for a metastatic 
% tumor growth model.
%
% This reproduces figure 3 from "Numerical solution of metastatic tumor 
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

fprintf('Welcome to MTGM demo #6\n');
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

tic
% No treatment 
fprintf('No treatment \n');
varargin_no_treat = {};
tic
[M_no_treat,N_no_treat,j_no_treat,C_no_treat] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_no_treat{:}); 
toc
fprintf(' Index j and condition number C with No treatment \n');
[j_no_treat, C_no_treat]


% Endostatine 1
fprintf('Endostatin 1 \n');
varargin_end_1 = {'treatment_type','end_1'};
tic
[M_end_1,N_end_1,j_end_1,C_end_1] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_1{:}); 
toc
fprintf(' Index j and condition number C with Endostatin 1 \n');
[j_end_1, C_end_1]

% Endostatin 2
fprintf('Endostatin 2 \n');
varargin_end_2 = {'treatment_type','end_2'};
tic
[M_end_2,N_end_2,j_end_2,C_end_2] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_2{:}); 
toc

fprintf(' Index j and condition number C with Endostatin 2 \n');
[j_end_2, C_end_2]

% Endostatine 3
fprintf('Endostatin 3 \n');
varargin_end_3 = {'treatment_type','end_3'};
tic
[M_end_3,N_end_3,j_end_3,C_end_3] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_3{:}); 
toc

fprintf(' Index j and condition number C with Endostatin 3 \n');
[j_end_3, C_end_3]

% change the directory to save the figures: 
cd Figures

% plot and save the results
fprintf('Plot metastatic mass \n');
figure
plot(t,M_no_treat,'b-*',t,M_end_1,'r-+',t,M_end_2,'k-o',...
    t,M_end_3,'m-^','linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
legend('Without treatment','Endostatin 1','Endostatin 2','Endostatin 3','Location','best');
axis([0 15 0 0.5])
saveas(gcf, 'Met_mass_demo6', 'fig');
saveas(gcf, 'Met_mass_demo6', 'epsc');

fprintf('Plot cumulative number of metastases \n');
figure
plot(t,N_no_treat,'b-*', t, N_end_1,'r-+',t,N_end_2,'k-o',t,...
    N_end_3,'m-^','linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Cumulative number')
legend('Without treatment','Endostatin 1','Endostatin 2','Endostatin 3','Location','best');
saveas(gcf, 'Cum_numb_demo6', 'fig');
saveas(gcf, 'Cum_numb_demo6', 'epsc');


