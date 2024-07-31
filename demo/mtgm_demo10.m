%--------------------------------------------------------------------------
% mtgm_demo10:
%
% This demo solves the VIE in order to compute the cumulative number of
% metastases, N(t), and the total metastatic mass M(t) for a metastatic 
% tumor growth model.
%
% This reproduces figure 8 from "A new MATLAB software for numerical
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

fprintf('Welcome to MTGM demo #10\n');
fprintf(['Compute the metastatic mass, M(t), and the cumulative number of ' ...
    'metastases, N(t) \n'])
fprintf('2D non-autonomous model without treatment \n')
addpath('..');
addpath('../utils');

% time in days
T = 360;
T0 = 15;
t = (1:T);

fprintf('No treatment for k and k-1 \n');
% change the initial condition of x0_p with respect to default values
varargin_no_treat = {'y0_p', [10^-6, 625]};
k = 6;
kind = 3;
type_OB = 0;

tic
[M_no_treat_0, N_no_treat_0] = ...
    VieSolve(kind,T0,t,2^k,type_OB,varargin_no_treat{:});
toc

tic
[M_no_treat, N_no_treat,j_no_treat,C_no_treat] = ...
    VieSolve(kind,T0,t,2^(k-1),type_OB,varargin_no_treat{:});
toc

% control of the condition number of the solved linear system
if C_no_treat > 1/eps
    sprintf('Warning: Results may be inaccurate. COND = %d',C_no_treat)
end

fprintf('Index j and condition number C with No treatment \n');
[j_no_treat, C_no_treat]
% compute the relative errors related to the approximate values of the
% Metastatic Mass
err_rel_M_no_treat = abs(M_no_treat_0-M_no_treat)./abs(M_no_treat_0);
% compute the relative errors related to the approximate values of the
% cumulative number of metasteses
err_rel_N_no_treat = abs(N_no_treat_0-N_no_treat)./abs(N_no_treat_0);


errM_no_treat= max(err_rel_M_no_treat);
errN_no_treat = max(err_rel_N_no_treat);
fprintf('Error for metastatic mass and for cumulative number of metastases \n');
[errM_no_treat, errN_no_treat]

% change the directory to save the figures: 
cd Figures

% plot and save the results
fprintf('Plot metastatic mass \n');
figure
plot(t,M_no_treat,'b-*','linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
axis([0 360 0 7*10^15]);
saveas(gcf, 'Met_mass_demo10', 'fig');
saveas(gcf, 'Met_mass_demo10', 'epsc');

fprintf('Plot cumulative number of metastases \n');
figure
plot(t,N_no_treat,'b-*','linewidth',2)
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Cumulative number')
axis([0 360 0 10*10^12]);
saveas(gcf, 'Cum_numb_demo10', 'fig');
saveas(gcf, 'Cum_numb_demo10', 'epsc');
