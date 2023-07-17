close all; clear; clc;

format long

tic

% Compute full equilibrium distributions for ten-step block-pulse model for a given rho value.
% Uses an Allee growth function

% Runs through each possible p-level equilibrium option and solves for the
% p-1 spatial-threshold r_i values by solving the nonlinear system of equations
% n(r_i) = iN/m = ih, i = q, ..., q+p-2 (where p is the number of levels in the
% equilibrium solution and q is the index of the smallest growth level 
% seen in the equilibrium).

%Parameters
L = 1;          %Domain length
N = 1;          %Upper population density limit
K = 1;          %Carrying capacity
alpha = 5;      %1/avg dispersal distance
rho = 2.8;      %Growth parameter
m = 10;         %Steps in block-pulse model
h = N/m;

%Set up spatial grid
dx_conv = 2^-10;
x = -L/2:dx_conv:L/2;

% Laplace PDF
lpdf = @(x) (1/2)*alpha*exp(-alpha*abs(x));
lCDF = @(x) (1/2).*(1 + sign(x).*(1 - exp(-alpha.*abs(x))));

%Define algorithm for use in solving nonlinear system of equations
options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg');

%% Initialize storage of p-rate equilibrium solutions

%One growth rate
n_one = cell(m, 1); 
%Two growth rates
n_two = cell(m-1, 2); 
%Three growth rates
n_three = cell(m-2, 2); 
%Four growth rates
n_four = cell(m-3, 2);
%Five growth rates, first solution
n_five = cell(m-4, 2);
%Six growth rates, first solution
n_six = cell(m-5, 2);
%Seven growth rates, first solution
n_seven = cell(m-6, 2);
%Seven growth rates, first solution
n_eight = cell(m-7, 2);
%Seven growth rates, first solution
n_nine = cell(m-8, 2);
%Seven growth rates, first solution
n_ten = cell(m-9, 2);

%% Block-pulse coefficients g_i
%For the Allee model
gcoefs = zeros(1, m);
%Solve block-pulse coefficients
for k = 1:m
    gcoefs(k) = integral(@(n) m*((((1 + rho^2)./K).*n.^2)./(1 + (rho/K)^2.*n.^2)), (k-1)*(1/m), k*(1/m));
end

%% One-level equilibria

%Calculate n(x)
for k = 1:m
    %Calculate n(0) and n(L/2)
    n0 = gcoefs(k)*(lCDF(L/2) - lCDF(L/2));
    nL2 = gcoefs(k)*(lCDF(L) - lCDF(0));
    %Sanity check: keep only solutions where n(0) < kN/m and n(L/2) > (k-1)N/m
    if n0 < k*N/m && nL2 > (k-1)*N/m
        n_one{k} = gcoefs(k)*(lCDF(x + L/2) - lCDF(x - L/2));
    else
        n_one{k} = NaN;
    end
end

%% Two-level equilibria

%Calculate n(x)
for j = 1:m-1
    %Construct n(r_j) = j*N/m function
    rjfunc = @(rj) gcoefs(j).*(lCDF(rj + L/2) - lCDF(rj - L/2)) + ...
        (gcoefs(j+1) - gcoefs(j)).*(lCDF(2.*rj) - lCDF(0)) - j*N/m;
    
    %Find r_j solution(s):
    r1 = fzero(rjfunc, 0.001);
    r2 = fzero(rjfunc, L/2);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 2, gcoefs(j:j+1), r);
    nL2 = @(r) pratenL2(L, 2, gcoefs(j:j+1), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > j*h && n0(r1) < (j+1)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            r1 > 0 && r1 < L/2
        n_two{j,1} = pratesol(L, 2, gcoefs(j:j+1), r1, x);
    else
        n_two{j,1} = NaN;
    end
    if n0(r2) > j*h && n0(r2) < (j+1)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            r2 > 0 && r2 < L/2
        n_two{j,2} = pratesol(L, 2, gcoefs(j:j+1), r2, x);
    else
        n_two{j,2} = NaN;
    end 
end

%% Three-level equilibria:

%Calculate n(x)
for j = 1:m-2
    %Construct n(r_j) = j*N/m and n(r_j+1) = (j+1)*N/m function pair
    rfunc = @(r) threeratesol(L, j, m, gcoefs(j), gcoefs(j+1), gcoefs(j+2), r);
    r_init1 = [L/2-0.01, L/2-0.01];
    r_init2 = [L/2-0.01, 0.01];
    
    %Find r_j solution(s):
    r1 = fsolve(rfunc, r_init1, options);
    r2 = fsolve(rfunc, r_init2, options);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 3, gcoefs(j:j+2), r);
    nL2 = @(r) pratenL2(L, 3, gcoefs(j:j+2), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > (j+1)*h && n0(r1) < (j+2)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            0 < r1(2) && r1(2) < r1(1) && r1(1) < L/2
        n_three{j,1} = pratesol(L, 3, gcoefs(j:j+2), r1, x);
    else
        n_three{j,1} = NaN;
    end
    if n0(r2) > (j+1)*h && n0(r2) < (j+2)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            0 < r2(2) && r2(2) < r2(1) && r2(1) < L/2
        n_three{j,2} = pratesol(L, 3, gcoefs(j:j+2), r2, x);
    else
        n_three{j,2} = NaN;
    end 
end

%% Four-level equilibria:

%Calculate n(x)
for j = 1:m-3
    %Construct n(r_j) = j*N/m, n(r_j+1) = (j+1)*N/m, n(r_j+2) = (j+2)*N/m function triple
    rfunc = @(r) fourratesol(L, j, m, gcoefs(j), gcoefs(j+1), gcoefs(j+2), gcoefs(j+3), r);
    r_init1 = [L/2-0.001, L/2-0.01, L/4];
    r_init2 = [L/4, L/4, 0.001];
    
    %Find r_j solution(s):
    r1 = fsolve(rfunc, r_init1, options);
    r2 = fsolve(rfunc, r_init2, options);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 4, gcoefs(j:j+3), r);
    nL2 = @(r) pratenL2(L, 4, gcoefs(j:j+3), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > (j+2)*h && n0(r1) < (j+3)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            0 < r1(3) && r1(3) < r1(2) && r1(2) < r1(1) && r1(1) < L/2
        n_four{j,1} = pratesol(L, 4, gcoefs(j:j+3), r1, x);
    else
        n_four{j,1} = NaN;
    end
    if n0(r2) > (j+2)*h && n0(r2) < (j+3)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            0 < r2(3) && r2(3) < r2(2) && r2(2) < r2(1) && r2(1) < L/2
        n_four{j,2} = pratesol(L, 4, gcoefs(j:j+3), r2, x);
    else
        n_four{j,2} = NaN;
    end 
end

%% Five-level equilibria:

%Calculate n(x)
for j = 1:m-4
    %Construct n(r_j) = j*N/m function quad
    rfunc = @(r) fiveratesol(L, j, m, gcoefs(j), gcoefs(j+1), gcoefs(j+2), gcoefs(j+3), gcoefs(j+4), r);
    r_init1 = [L/2 - 0.01, L/4, L/4, L/4];
    r_init2 = [L/4, L/4, L/4, 0.001];
    
    %Find r_j solution(s):
    r1 = fsolve(rfunc, r_init1, options);
    r2 = fsolve(rfunc, r_init2, options);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 5, gcoefs(j:j+4), r);
    nL2 = @(r) pratenL2(L, 5, gcoefs(j:j+4), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > (j+3)*h && n0(r1) < (j+4)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            0 < r1(4) && r1(4) < r1(3) && r1(3) < r1(2) && r1(2) < r1(1) && r1(1) < L/2
        n_five{j,1} = pratesol(L, 5, gcoefs(j:j+4), r1, x);
    else
        n_five{j,1} = NaN;
    end
    if n0(r2) > (j+3)*h && n0(r2) < (j+4)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            0 < r2(4) && r2(4) < r2(3) && r2(3) < r2(2) && r2(2) < r2(1) && r2(1) < L/2
        n_five{j,2} = pratesol(L, 5, gcoefs(j:j+4), r2, x);
    else
        n_five{j,2} = NaN;
    end 
end

%% Six-level equilibria:

%Calculate n(x)
for j = 1:m-5
    %Construct n(r_j) = j*N/m function quint
    rfunc = @(r) sixratesol(L, j, m, gcoefs(j), gcoefs(j+1), gcoefs(j+2), gcoefs(j+3), gcoefs(j+4), gcoefs(j+5), r);
    r_init1 = [L/2 - 0.01, L/4, L/4, L/4, L/4];
    r_init2 = [L/4, L/4, L/4, L/4, 0.001];
    
    %Find r_j solution(s):
    r1 = fsolve(rfunc, r_init1, options);
    r2 = fsolve(rfunc, r_init2, options);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 6, gcoefs(j:j+5), r);
    nL2 = @(r) pratenL2(L, 6, gcoefs(j:j+5), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > (j+4)*h && n0(r1) < (j+5)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            0 < r1(5) && r1(5) < r1(4) && r1(4) < r1(3) && r1(3) < r1(2) && r1(2) < r1(1) && r1(1) < L/2
        n_six{j,1} = pratesol(L, 6, gcoefs(j:j+5), r1, x);
    else
        n_six{j,1} = NaN;
    end
    if n0(r2) > (j+4)*h && n0(r2) < (j+5)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            0 < r2(5) && r2(5) < r2(4) && r2(4) < r2(3) && r2(3) < r2(2) && r2(2) < r2(1) && r2(1) < L/2
        n_six{j,2} = pratesol(L, 6, gcoefs(j:j+5), r2, x);
    else
        n_six{j,2} = NaN;
    end 
end

%% Seven-level equilibria:

%Calculate n(x)
for j = 1:m-6
    %Construct n(r_j) = j*N/m function quint
    rfunc = @(r) sevenratesol(L, j, m, gcoefs(j), gcoefs(j+1), gcoefs(j+2), gcoefs(j+3), ...
        gcoefs(j+4), gcoefs(j+5), gcoefs(j+6), r);
    r_init1 = [L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01];
    r_init2 = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001];
    
    %Find r_j solution(s):
    r1 = fsolve(rfunc, r_init1, options);
    r2 = fsolve(rfunc, r_init2, options);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 7, gcoefs(j:j+6), r);
    nL2 = @(r) pratenL2(L, 7, gcoefs(j:j+6), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > (j+5)*h && n0(r1) < (j+6)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            0 < r1(6) && r1(6) < r1(5) && r1(5) < r1(4) && r1(4) < r1(3) && ...
                r1(3) < r1(2) && r1(2) < r1(1) && r1(1) < L/2
        n_seven{j,1} = pratesol(L, 7, gcoefs(j:j+6), r1, x);
    else
        n_seven{j,1} = NaN;
    end
    if n0(r2) > (j+5)*h && n0(r2) < (j+6)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            0 < r2(6) && r2(6) < r2(5) && r2(5) < r2(4) && r2(4) < r2(3) && ...
                r2(3) < r2(2) && r2(2) < r2(1) && r2(1) < L/2
        n_seven{j,2} = pratesol(L, 7, gcoefs(j:j+6), r2, x);
    else
        n_seven{j,2} = NaN;
    end 
end

%% Eight-level equilibria:

%Calculate n(x)
for j = 1:m-7
    %Construct n(r_j) = j*N/m function sept
    rfunc = @(r) eightratesol(L, j, m, gcoefs(j), gcoefs(j+1), gcoefs(j+2), gcoefs(j+3), ...
        gcoefs(j+4), gcoefs(j+5), gcoefs(j+6), gcoefs(j+7), r);
    r_init1 = [L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01];
    r_init2 = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001];
    
    %Find r_j solution(s):
    r1 = fsolve(rfunc, r_init1, options);
    r2 = fsolve(rfunc, r_init2, options);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 8, gcoefs(j:j+7), r);
    nL2 = @(r) pratenL2(L, 8, gcoefs(j:j+7), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > (j+6)*h && n0(r1) < (j+7)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            0 < r1(7) && r1(7) < r1(6) && r1(6) < r1(5) && r1(5) < r1(4) && r1(4) < r1(3) && ...
                r1(3) < r1(2) && r1(2) < r1(1) && r1(1) < L/2
        n_eight{j,1} = pratesol(L, 8, gcoefs(j:j+7), r1, x);
    else
        n_eight{j,1} = NaN;
    end
    if n0(r2) > (j+6)*h && n0(r2) < (j+7)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            0 < r2(7) && r2(7) < r2(6) && r2(6) < r2(5) && r2(5) < r2(4) && r2(4) < r2(3) && ...
                r2(3) < r2(2) && r2(2) < r2(1) && r2(1) < L/2
        n_eight{j,2} = pratesol(L, 8, gcoefs(j:j+7), r2, x);
    else
        n_eight{j,2} = NaN;
    end 
end

%% Nine-level equilibria:

%Calculate n(x)
for j = 1:m-8
    %Construct n(r_j) = j*N/m function oct
    rfunc = @(r) nineratesol(L, j, m, gcoefs(j), gcoefs(j+1), gcoefs(j+2), gcoefs(j+3), ...
        gcoefs(j+4), gcoefs(j+5), gcoefs(j+6), gcoefs(j+7), gcoefs(j+8), r);
    r_init1 = [L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01];
    r_init2 = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001];
    
    %Find r_j solution(s):
    r1 = fsolve(rfunc, r_init1, options);
    r2 = fsolve(rfunc, r_init2, options);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 9, gcoefs(j:j+8), r);
    nL2 = @(r) pratenL2(L, 9, gcoefs(j:j+8), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > (j+7)*h && n0(r1) < (j+8)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            0 < r1(8) && r1(8) < r1(7) && r1(7) < r1(6) && r1(6) < r1(5) && r1(5) < r1(4) && r1(4) < r1(3) && ...
                r1(3) < r1(2) && r1(2) < r1(1) && r1(1) < L/2
        n_nine{j,1} = pratesol(L, 9, gcoefs(j:j+8), r1, x);
    else
        n_nine{j,1} = NaN;
    end
    if n0(r2) > (j+7)*h && n0(r2) < (j+8)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            0 < r2(8) && r2(8) < r2(7) && r2(7) < r2(6) && r2(6) < r2(5) && r2(5) < r2(4) && r2(4) < r2(3) && ...
                r2(3) < r2(2) && r2(2) < r2(1) && r2(1) < L/2
        n_nine{j,2} = pratesol(L, 9, gcoefs(j:j+8), r2, x);
    else
        n_nine{j,2} = NaN;
    end 
end

%% Ten-level equilibria:

%Calculate n(x)
for j = 1:m-9
    %Construct n(r_j) = j*N/m function nona
    rfunc = @(r) tenratesol(L, j, m, gcoefs(j), gcoefs(j+1), gcoefs(j+2), gcoefs(j+3), ...
        gcoefs(j+4), gcoefs(j+5), gcoefs(j+6), gcoefs(j+7), gcoefs(j+8), gcoefs(j+9), r);
    r_init1 = [L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01, L/2 - 0.01];
    r_init2 = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001];
    
    %Find r_j solution(s):
    r1 = fsolve(rfunc, r_init1, options);
    r2 = fsolve(rfunc, r_init2, options);
    
    %Calculate n0 and nL2
    n0 = @(r) praten0(L, 10, gcoefs(j:j+9), r);
    nL2 = @(r) pratenL2(L, 10, gcoefs(j:j+9), r);
    
    %Keep solutions that satisfy a sanity check
    if n0(r1) > (j+8)*h && n0(r1) < (j+9)*h && ...
            nL2(r1) > (j-1)*h && nL2(r1) < j*h && ...
            0 < r1(9) && r1(9) < r1(8) && r1(8) < r1(7) && r1(7) < r1(6) && r1(6) < r1(5) && ...
            r1(5) < r1(4) && r1(4) < r1(3) && r1(3) < r1(2) && r1(2) < r1(1) && r1(1) < L/2
        n_ten{j,1} = pratesol(L, 10, gcoefs(j:j+9), r1, x);
    else
        n_ten{j,1} = NaN;
    end
    if n0(r2) > (j+8)*h && n0(r2) < (j+9)*h && ...
            nL2(r2) > (j-1)*h && nL2(r2) < j*h && ...
            0 < r2(9) && r2(9) < r2(8) && r2(8) < r2(7) && r2(7) < r2(6) && r2(6) < r2(5) && ...
            r2(5) < r2(4) && r2(4) < r2(3) && r2(3) < r2(2) && r2(2) < r2(1) && r2(1) < L/2
        n_ten{j,2} = pratesol(L, 10, gcoefs(j:j+9), r2, x);
    else
        n_ten{j,2} = NaN;
    end 
end

%% Plot equilibrium distributions

%rho = 2.8;
f1 = figure(1);
set(f1, 'Position', [10 400 850 475]);
ax(1) = subplot(1,2,1);
ax(1).Position = [0.1 0.15 0.5 0.8];
hold on
set(ax(1), 'Fontsize', 20, 'xlim', [-L/2, L/2], ...
    'ytick', [0 2*N/10 4*N/10 6*N/10 8*N/10 N], 'yticklabels', {0 2*N/10 4*N/10 6*N/10 8*N/10 N}, ...
    'TickDir', 'out', 'TickLength', [0.02 0.02]);
hold on
%Manually plots the only equilibria found for this rho value
plot(x, n_one{1}, 'color', [234/255 123/255 123/255], 'linewidth', 2);
plot(x, n_three{1,2}, 'color', [70/255 70/255 236/255], 'linewidth', 2);
plot(x, n_five{5,1}, 'color', [230/255 25/255 28/255], 'linewidth', 2);
%This option plots all equilibria by cycling through each cell rather than
%manually sorting through equilibria.
% %One-level equilibria
% cellfun(@(y) plot(x, y, 'color', [234/255 123/255 123/255], 'linewidth', 2), n_one)
% %Two-level equilibria
% cellfun(@(y) plot(x, y, 'color', [25/255 205/255 225/255], 'linewidth', 2), n_two)
% %Three-level equilibria
% cellfun(@(y) plot(x, y, 'color', [70/255 70/255 236/255], 'linewidth', 2), n_three)
% %Four-level equilibria
% cellfun(@(y) plot(x, y, 'color', [243/255 184/255 97/255], 'linewidth', 2), n_four)
% %Five-level equilibria
% cellfun(@(y) plot(x, y, 'color', [230/255 25/255 28/255], 'linewidth', 2), n_five)
% %Six-level equilibria
% cellfun(@(y) plot(x, y, 'color', [11/255 218/255 81/255], 'linewidth', 2), n_six)
% %Seven-level equilibria
% cellfun(@(y) plot(x, y, 'color', [153/255, 50/255, 204/255], 'linewidth', 2), n_seven)
% %Eight-level equilibria
% cellfun(@(y) plot(x, y, 'color', [0/255 0/255 0/255], 'linewidth', 2), n_eight)
% %Nine-level equilibria
% cellfun(@(y) plot(x, y, 'color', [207/255 0/255 222/255], 'linewidth', 2), n_nine)
% %Ten-level equilibria
% cellfun(@(y) plot(x, y, 'color', [12/255 89/255 0/255], 'linewidth', 2), n_ten)

%Add density thresholds:
for i = 1:m-1
    yline(i*N/m, 'k--');
end
xlim([-L/2, L/2]);
xlabel('x (space)'); ylabel('n(x)');

%Legend
Lega(1) = plot(nan, nan, 'color', [234/255 123/255 123/255], 'linestyle', '-', 'linewidth', 3);
Lega(2) = plot(nan, nan, 'color', [25/255 205/255 225/255], 'linestyle', '-', 'linewidth', 3);
Lega(3) = plot(nan, nan, 'color', [70/255 70/255 236/255], 'linestyle', '-', 'linewidth', 3);
Lega(4) = plot(nan, nan, 'color', [243/255 184/255 97/255], 'linestyle', '-', 'linewidth', 3);
Lega(5) = plot(nan, nan, 'color', [230/255 25/255 28/255], 'linestyle', '-', 'linewidth', 3);
Lega(6) = plot(nan, nan, 'color', [11/255 218/255 81/255], 'linestyle', '-', 'linewidth', 3);
Lega(7) = plot(nan, nan, 'color', [153/255, 50/255, 204/255], 'linestyle', '-', 'linewidth', 3);
Lega(8) = plot(nan, nan, 'color', [0/255 0/255 0/255], 'linestyle', '-', 'linewidth', 3);
Lega(9) = plot(nan, nan, 'color', [207/255 0/255 222/255], 'linestyle', '-', 'linewidth', 3);
Lega(10) = plot(nan, nan, 'color', [12/255 89/255 0/255], 'linestyle', '-', 'linewidth', 3);
legend(Lega, {'One-level solution', 'Two-level solution', 'Three-level solution', 'Four-level solution',...
    'Five-level solution', 'Six-level solution', 'Seven-level solution', 'Eight-level solution', ...
    'Nine-level solution', 'Ten-level solution',}, ...
    'Position', [0.71 0.4 0.2 0.3])
