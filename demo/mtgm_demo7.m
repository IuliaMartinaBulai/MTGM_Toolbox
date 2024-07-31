%--------------------------------------------------------------------------
% mtgm_demo7:
%
% This demo solves the VIE in order to compare the effect of five different
% doses of Endostatin on the cumulative number of metastases, N(t), and on 
% the total metastatic mass, M(t).
%
% This reproduces figure 4 from "Numerical solution of metastatic tumor 
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

fprintf('Welcome to MTGM demo #7\n');
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
fprintf('No treatment \n');
varargin_no_treat = {};
tic
[M_no_treat,N_no_treat,j_no_treat,C_no_treat] = ...
   VieSolve(kind,T0,t,2^k,type_OB,varargin_no_treat{:}); 
toc
fprintf(' Index j and condition number C with No treatment \n');
[j_no_treat, C_no_treat]

% Endostatine 20 mg
fprintf('Endostatine 20 mg \n');
varargin_end_1 = {'treatment_type','end_1'};
tic
[M_end_1,N_end_1,j_end_1,C_end_1] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_1{:}); 
toc
fprintf(' Index j and condition number C with Endostatine 20 mg \n');
[j_end_1, C_end_1]

% Endostatin 5 mg
fprintf('Endostatine 5 mg \n');
varargin_end_2 = {'treatment_type','end_4'};
tic
[M_end_2,N_end_2,j_end_2,C_end_2] =...
   VieSolve(kind,T0,t,2^k,type_OB,varargin_end_2{:}); 
toc
fprintf(' Index j and condition number C with Endostatine 5 mg \n');
[j_end_2, C_end_2]

% Endostatine 10 mg
fprintf('Endostatine 10 mg \n');
varargin_end_3 = {'treatment_type','end_5'};
tic
[M_end_3,N_end_3,j_end_3,C_end_3] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_3{:}); 
toc
fprintf(' Index j and condition number C with Endostatine 10 mg \n');
[j_end_3, C_end_3]

% Endostatine 30 mg
fprintf('Endostatine 30 mg \n');
varargin_end_4 = {'treatment_type','end_6'};
tic
[M_end_4,N_end_4,j_end_4,C_end_4] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_4{:}); 
toc
fprintf(' Index j and condition number C with Endostatine 30 mg \n');
[j_end_4, C_end_4]

% Endostatine 40 mg
fprintf('Endostatine 40 mg \n');
varargin_end_5 = {'treatment_type','end_7'};
tic
[M_end_5,N_end_5,j_end_5,C_end_5] =...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_5{:}); 
toc
fprintf(' Index j and condition number C with Endostatine 40 mg \n');
[j_end_5, C_end_5]

% change the directory to save the figures: 
cd Figures

% plot and save the results
fprintf('Plot metastatic mass \n');
figure
plot(t,M_no_treat,'b-*',t,M_end_2,'r-+',t,M_end_3,'k-o',...
    t,M_end_1,'m-^',t,M_end_4,'g-s',t,M_end_5,'r-*', 'linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
legend('Without treatment','Dose = 5mg','Dose = 10mg','Dose = 20mg',...
    'Dose = 30mg','Dose = 40mg','Location','best');
axis([0 15 0 1.5]);
saveas(gcf, 'Met_mass_demo7', 'fig');
saveas(gcf, 'Met_mass_demo7', 'epsc');


fprintf('Plot cumulative number of metastases \n');
figure
plot(t,N_no_treat,'b-*',t,N_end_2,'r-+',t,N_end_3,'k-o',...
    t,N_end_1,'m-^',t,N_end_4,'g-s',t,N_end_5,'r-*', 'linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Cumulative number')
legend('Without treatment','Dose = 5mg','Dose = 10mg','Dose = 20mg',...
    'Dose = 30mg','Dose = 40mg','Location','best');
saveas(gcf, 'Cum_numb_demo7', 'fig');
saveas(gcf, 'Cum_numb_demo7', 'epsc');

