% This script aims to compare the diffrent choices that can be made for a
% ResQPASS implementation using quadprog. These results are also compare
% with the performance of the 3-th version of ResQPASS from  
% https://github.com/AppliedMathUAntwerpen/ResQPASS/tree/main/Matlab

% This corresponds to Example 4.3.1 (Figure 4.6, p. 28) in the accompanying paper.
clear all; close all; clc;
rng(42);

% Determine whether the results chould be saved
saveFigure = false;
%% %%%%%%%%%%%%%%%%%%%%% Initialize test problem %%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomly generate the matrix A and vector b
m = 300;
n = 100;
A = 0.15 * rand(m, n);
b = randn(m, 1);

% Set the lower and upper bounds
lb = zeros(n, 1) - 1;
ub = zeros(n, 1) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%% Solve the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WarmTimes = zeros(100, 1);
NoWarmTimes = zeros(100, 1);
intPoint = zeros(100, 1);
original_NoWarmStart = zeros(100, 1);
original_WarmStart = zeros(100, 1);

% Solve the problem using diffrent techniques
for k = 1:100
    tic
    [~, history_warmStart] = ResQPASS_quadprog(A, b, lb, ub, "warmStart", true);        % With warm starting
    WarmTimes(k) = toc;
    
    tic
    [~, history_noWarmStart] = ResQPASS_quadprog(A, b, lb, ub, "warmStart", false);     % Without warm starting
    NoWarmTimes(k) = toc;
    
    tic
    [~, ~] = ResQPASS_quadprog(A, b, lb, ub, "algorithm", "interior-point-convex");     % Interior point
    intPoint(k) = toc;

    tic
    [~, history_orig_NWS] = ResQPASS_altered(A, b, lb, ub, "doWarmStart", false);       % Original algorithm, no warm starting
    original_NoWarmStart(k) = toc;

    tic 
    [~, history_orig_WS] = ResQPASS_altered(A, b, lb, ub, ...
        "doWarmStart", true, "maxInnerIt", max(history_orig_NWS.innerIterations));      % Original algorithm, with warm starting
    original_WarmStart(k) = toc;
end
%% Plot the amount of inneriterations
figure;

% Define colors for consistency.
colors = lines(5);

%%%%%%%%%%%%%%%%%%%%%%%
% Time each test took %
%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 2, 1);
hold on;

plot(NoWarmTimes,            'Color', colors(1,:))
plot(WarmTimes,              'Color', colors(2,:))
plot(intPoint,               'Color', colors(3,:))
plot(original_NoWarmStart,   'Color', colors(4,:))
plot(original_WarmStart,     'Color', colors(5,:))

xlabel("Test nr")
ylabel("Tijd tot convergentie (s)")
legend("Active set", "Active set, warme start", ...
    "Interior point", "ResQPASS", ...
    "ResQPASS, warme start", "Interpreter", "latex") 
title("Convergentietijd implementaties")
ylim([0, max(intPoint) + 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of inner iterations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 2, 2);
hold on;

plot(history_noWarmStart.innerIterations, 'Color', colors(1,:))
plot(history_warmStart.innerIterations,   'Color', colors(2,:))
plot(history_orig_NWS.innerIterations,    'Color', colors(4,:))
plot(history_orig_WS.innerIterations,     'Color', colors(5,:))

active_bounds = nnz(history_warmStart.lambda(:, end) > 0) + ...
                nnz(history_warmStart.mu(:, end) > 0);
xline(active_bounds, "g--")

xlabel('Iteratie externe loop')
ylabel('Iteraties interne loop')
legend("Active set", 'Active set, warme start', ...
    "ResQPASS", ...
    "ResQPASS, warme start", "$\vert \mathcal{W} \vert$", ...
    "Interpreter", "latex")
title("Effect van warme start")

if saveFigure, setFigParameters("warm_starting"); end