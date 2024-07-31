%----------------------------------------------------------------------------------------
% File: Gbeta.m
%
% Goal: compute the right-hand side of the VIE
%
% Use: [gtbeta] = Gbeta(kind,t,w1,Y,type_OB,varargin)
%
% Input: kind - see VieSolve.m file
%        t - row array of the evaluation points
%        w1 - weights of the M-point Gauss-Laguerre rule with with M = 2048
%        Y - element-wise product between the array of the evaluation points and the
%            array exp(-x1) with x1 array of the zeros of the M-th Laguerre polynomial
%        type_OB - see VieSolve.m
%        varargin - see VieSolve.m
%
% Output: gtbeta - array of the right-hand side of the VIE at the evaluation points t
%
% Recalls: vie.m, tumorGrowthFunctions.m, tumorGrowthFunctionsOdeSolve.m, ode23t.m,
%          odefun1Dm .m, odefun1Dp.m tumorTreatmentParam.m, odefun2Dm.m, odefun2Dp.m
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
function [gtbeta] = Gbeta(kind,t,w1,Y,type_OB,varargin)

% kind = 0: Volterra Integral Equation
if kind == 0
    [g] = vie();
    g2 = g{2};        
    gtbeta = g2(t);
% kind = 1: Cumulative number of metastases and total metastatic mass for
%           the 1D metastatic tumor growth model
elseif kind == 1
    [g,gamma,~,~] = tumorGrowthFunctions(varargin{:});

    if type_OB == 1
        h1 = g{1};
        m = length(t);
        gtbeta = zeros(1,m);
        for i = 1:m
            y = t(i)+gamma;
            r = 1;
            g = w1(r)*y*h1(Y(r,i),t(i));
            while (r< length(w1)) && (abs(g)>0.5e-25)
                gtbeta(i) = gtbeta(i)+g;
                r = r+1;
                g = w1(r)*y*h1(Y(r,i),t(i));
            end
        end
    elseif type_OB == 2
        h2 = g{2};
        m = length(t);
        gtbeta = zeros(1,m);
        for i = 1:m
            y = t(i)+gamma;
            r = 1;
            g = w1(r)*y*h2(Y(r,i),t(i));
            while (r< length(w1)) && (abs(g)>0.5e-25)
                gtbeta(i) = gtbeta(i)+g;
                r = r+1;
                g = w1(r)*y*h2(Y(r,i),t(i));
            end
        end
    end
    % kind = 2: Cumulative number of metastases and total metastatic mass
    %           for the 1D metastatic tumor growth model solving the 1D ODE
elseif kind == 2
    [g,y0_m,y0_p,param,options] = tumorGrowthFunctionsOdeSolve(varargin{:});
    % type_OB = 1: compute the total metastatic mass
    if type_OB == 1
        be_p = g{1};
        m = length(t);
        gtbeta = zeros(1,m);
        parfor i = 1:m
            y = t(i);
            r = 1;
            [~,sol1] = ode45(@(tt,yy)odefun1Dm (tt,yy,param),[0 t(i) - Y(r,i)],y0_m,options);
            [~,sol2] = ode45(@(tt,yy)odefun1Dp(tt,yy,param),[0 Y(r,i)],y0_p,options);
            g1 = w1(r)*y*sol1(end,1)*be_p(sol2(end,1));
            while (r< length(w1)) && (abs(g1)>0.5e-25)
                gtbeta(i) = gtbeta(i)+g1;
                r = r+1;
                [~,sol1] = ode45(@(tt,yy)odefun1Dm (tt,yy,param),[0 t(i)-Y(r,i)],y0_m,options);
                [~,sol2] = ode45(@(tt,yy)odefun1Dp(tt,yy,param),[0 Y(r,i)],y0_p,options);
                g1 = w1(r)*y*sol1(end,1)*be_p(sol2(end,1));
            end
        end
        % type_OB = 2: compute the cumulative number of metastases
    elseif type_OB == 2
        be_p = g{1};
        m = length(t);
        gtbeta = zeros(1,m);
        parfor i = 1:m
            y = t(i);
            r = 1;
            [~,sol2]=ode45(@(tt,yy)odefun1Dp(tt,yy, param),[0 Y(r,i)],y0_p,options); %Y(r,i)
            g = w1(r)*y*be_p(sol2(end,1));
            while (r< length(w1)) && (abs(g)>0.5e-25)
                gtbeta(i) = gtbeta(i)+g;
                r = r+1;
                [~,sol2]=ode45(@(tt,yy)odefun1Dp(tt,yy,param),[0 Y(r,i)],y0_p,options);
                g = w1(r)*y*be_p(sol2(end,1));
            end
        end
    end
    % kind = 3: Cumulative number of metastases and total metastatic mass for
    %           the 2D metastatic tumor growth model with treatment solving
    %           the 2D ODE non-autonomous system
    elseif kind == 3 
        [g,y0_p,y0_m,param,~,options] = tumorTreatmentParam(varargin{:});
        be_p = g{1};
        be_m = g{2};
        m = length(t);
        gtbeta = zeros(1,m);
        parfor i = 1:m
            y = t(i);
            if y ~= 0
                r = 1;
                [~,sol1] = ode23t(@(tt,yy)odefun2Dm(tt, yy, param),[Y(r,i) t(i)],y0_m,options);
                [~,sol2] = ode23t(@(tt,yy)odefun2Dp(tt, yy, param),[0 Y(r,i)],y0_p,options);
                g1 = w1(r)*y*be_m(sol1(end,1))*be_p(sol2(end,1));
                while (r< length(w1)) && (abs(g1)>0.5e-25)
                    gtbeta(i) = gtbeta(i)+g1;
                    r = r+1;
                    [~,sol1]= ode23t(@(tt,yy)odefun2Dm(tt, yy, param),[Y(r,i) t(i)],y0_m,options);
                    [~,sol2]= ode23t(@(tt,yy)odefun2Dp(tt, yy, param),[0 Y(r,i)],y0_p,options);
                    g1 = w1(r)*y*be_m(sol1(end,1))*be_p(sol2(end,1));
                end
            end
        end
        % kind = 4: Cumulative number of metastases and total metastatic mass for
        %           the 2D metastatic tumor growth model solving the 2D ODE
        %           autonomous system.
    elseif kind == 4
        [g,y0_p,y0_m,param,~,options] = tumorTreatmentParam(varargin{:});

        if type_OB == 1
            be_p = g{1};
            m = length(t);
            gtbeta = zeros(1,m);
            parfor i = 1:m
                y = t(i);
                r = 1;
                [~,sol1] = ode23t(@(tt,yy)odefun2Dm(tt,yy,param),[0 t(i) - Y(r,i)],y0_m,options);
                [~,sol2] = ode23t(@(tt,yy)odefun2Dp(tt,yy,param),[0 Y(r,i)],y0_p,options);
                g1 = w1(r)*y*sol1(end,1)*be_p(sol2(end,1));
                while (r< length(w1)) && (abs(g1)>0.5e-25)
                    gtbeta(i) = gtbeta(i)+g1;
                    r = r+1;
                    [~,sol1] = ode23t(@(tt,yy)odefun2Dm(tt,yy,param),[0 t(i) - Y(r,i)],y0_m,options); 
                    [~,sol2] = ode23t(@(tt,yy)odefun2Dp(tt,yy,param),[0 Y(r,i)],y0_p,options);
                    g1 = w1(r)*y*sol1(end,1)*be_p(sol2(end,1));
                end
            end
        elseif type_OB == 2
            be_p = g{1};
            m = length(t);
            gtbeta = zeros(1,m);
            parfor i = 1:m
                y = t(i);
                r = 1;
                [~,sol2]=ode23t(@(tt,yy)odefun2Dp(tt,yy, param),[0 Y(r,i)],y0_p,options); %Y(r,i)
                g = w1(r)*y*be_p(sol2(end,1));
                while (r< length(w1)) && (abs(g)>0.5e-25)
                    gtbeta(i) = gtbeta(i)+g;
                    r = r+1;
                    [~,sol2]=ode23t(@(tt,yy)odefun2Dp(tt,yy,param),[0 Y(r,i)],y0_p,options);
                    g = w1(r)*y*be_p(sol2(end,1));
                end
            end
        end
end
end
