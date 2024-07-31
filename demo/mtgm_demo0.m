%--------------------------------------------------------------------------
% mtgm_demo0:
%
% This demo computes the weighted solution of a VIE
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

fprintf('Welcome to VIE demo #0\n');
addpath('..');
addpath('../utils');
fprintf('Solve the VIE  \n');

T = 30;
T0 = 30;
t =(0:T);
kind = 0;
type_OB = 0;

% % exact solution
f = @(t) t;
% the weight u 
u = @(y) y.^(1/4).*exp(-y/2);
% row array of the evaluation points
varargin = {};
k = 9;
for i = 2:k-1
    tic
    [M,~,j(i-1),C(i-1)] = VieSolve(kind,T0,t,2^k,type_OB,varargin{:});
    toc
    err(i-1) = max(abs((f(t).*u(t))'-M));
end
fprintf('Index j, condition number and error\n')
[j' C' err']


% change the directory to save the figures: 
cd Figures

% plot and save the results
fprintf('Plot the Nystrom interpolating function \n');
figure
plot(t,M, 'linewidth',2)
set(gca,'fontsize',16)
title('Weighted Nystrom interpolating function')
xlabel('t')
ylabel('u(t)f_m(t)')
saveas(gcf, 'vie_demo0', 'fig');
saveas(gcf, 'vie_demo0', 'epsc');

