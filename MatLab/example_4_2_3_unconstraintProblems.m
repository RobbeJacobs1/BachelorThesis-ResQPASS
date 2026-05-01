% In this test, the residuals of ResQPASS will be compared to those of LSQR
% in an unconstrained case.

% This corresponds to Example 4.2.3 (Figure 4.5, p. 26) in the accompanying paper.
clear all; close all; clc;
rng(42);

% Determine whether the results should be saved
saveFigure = true;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Initialise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate matrices A and b
m = 300;
n = 100;
A = 0.15 * rand(m, n);
b = randn(m, 1);

% Disable the boundaries
lb = -Inf * ones(n, 1);
ub = Inf * ones(n, 1);

%% Calculate solutions
% Solve the problem using ResQPASS
[ResQPASS_sol, ResQPASS_history] = ResQPASS_quadprog(A, b, lb, ub);
% Solve the problem using LSQR
[lsqr_sol, lsqr_history] = lsqr_self(A, b, zeros(n, 1), 1000, 1e-8);
% Solve the problem using ResQPASS as introduced in https://arxiv.org/pdf/2302.13616 
[RSol, RHistory] = ResQPASS_altered(A, b, lb, ub);

%% %%%%%%%%%%%%%%%% Determine related parameters to convergence %%%%%%%%%%%
figure;
%%%%%%%%%%%%%
% Residuals %
%%%%%%%%%%%%%
subplot(1, 2, 1);
hold on
yscale("log");

% Plot the residuals
semilogy(ResQPASS_history.residuals/ResQPASS_history.residuals(1), "bo-");    % ResQPASS
semilogy(RHistory.residuals/RHistory.residuals(1))
semilogy(lsqr_history.residuals/lsqr_history.residuals(1), "rx-"); % LSQR

% Add aditional plot settings like legend, title, labels
legend("ResQPASS quadprog", "ResQPASS", "LSQR");                   % Legend
title("Residuen van LSQR en ResQPASS");                 % Title
xlabel('Iteratie $k$', "Interpreter", "Latex");         % x-label
ylabel("$\Vert r_k \Vert$", "Interpreter", "latex");    % y-label

%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the solutions %
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 2, 2);
hold on;

% Plot the solutions
plot(ResQPASS_sol, "bx-");      % ResQPASS
plot(lsqr_sol, "r--o");         % LSQR

% Add extra plot settings (title, legend, labels)
title("Beide oplossingen");                 % Title
legend("ResQPASS", "LSQR");                 % Legend
xlabel("$i$", "Interpreter", "Latex");      % x-label
ylabel("$x_i$", "Interpreter", "latex")     % y-label

% Save the results
if saveFigure setFigParameters("lsqr_vs_ResQPASS"); end