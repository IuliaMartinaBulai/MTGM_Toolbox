%--------------------------------------------------------------------------
% mtgm_demo9:
%
% This demo solves the VIE in order to compare the effect of three
% different drugs on the cumulative number of metastases, N(t), and on the
% total metastatic mass, M(t).
%
% This reproduces figure 7 from "A new MATLAB software for numerical
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

fprintf('Welcome to MTGM demo #9\n');
fprintf(['Compute the metastatic mass, M(t), and the cumulative number of ' ...
    'metastases, N(t) \n'])
fprintf('2D non-autonomous model with treatment \n')
addpath('..');
addpath('../utils');

% time in days

T = 360;
T0 = 15;
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


% Endostatine 1
fprintf('Endostatin 1 \n');
varargin_end_1 = {'treatment_type','end_1'};
tic
[M_end_1,N_end_1,j_end_1,C_end_1] = ...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_end_1{:});
toc
fprintf(' Index j and condition number C with Endostatin 1 \n');
[j_end_1, C_end_1]


% TNP_470
fprintf('TNP-470 \n');
varargin_TNP_470 = {'treatment_type','TNP-470'};
tic
[M_TNP_470,N_TNP_470,j_TNP_470,C_TNP_470]...
    = VieSolve(kind,T0,t,2^k,type_OB,varargin_TNP_470{:});
toc
fprintf(' Index j and condition number C with TNP-470 \n');
[j_TNP_470, C_TNP_470]


% Angiostatin
fprintf('angiostatin \n');
varargin_angiostatin = {'treatment_type','angiostatin'};
tic
[M_angiostatin,N_angiostatin,j_angiostatin,C_angiostatin] = ...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_angiostatin{:});
toc
fprintf(' Index j and condition number C with Angiostatin \n');
[j_angiostatin, C_angiostatin]

% change the directory to save the figures: 
cd Figures

% plot and save the results
fprintf('Plot metastatic mass \n');
figure
plot(t,M_TNP_470,'r-+',t,M_end_1,'k-o', t,M_angiostatin,'m-^',...
    t,M_no_treat,'b-*', 'linewidth',2);
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
legend('TPN-470','Endostatin','Angiostatin','Without treatment','Location','best');
axis([0 360 0 6*10^16]);
saveas(gcf, 'Met_mass_demo9', 'fig');
saveas(gcf, 'Met_mass_demo9', 'epsc');

fprintf('Plot cumulative number of metastases \n');
figure
plot(t, N_TNP_470,'r-+',t,N_end_1,'k-o', t, N_angiostatin,'m-^',...
    t,N_no_treat,'b-*','linewidth',2);
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Cumulative number')
legend('TPN-470','Endostatin','Angiostatin','Without treatment','Location','best');
axis([0 360 0 3.5*10^13]);
saveas(gcf, 'Cum_numb_demo9', 'fig');
saveas(gcf, 'Cum_numb_demo9', 'epsc');
