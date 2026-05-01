% In this script, the ResQPASS-algorithm is applied to a first randomly
% generated test problem. This problem is solved using the MatLab-function
% quadprog, the ResQPASS algorithm and in the unconstraint case. This is
% used to compare the results of the diff rent algorithms and examine
% convergence behaviour.

% This corresponds to Example 4.1.1 (Figures 4.1 and 4.2, p. 23) in the accompanying paper.
clear all; close all; clc;
rng(42);

% Determine whether the results should be saved.
saveFigure = false;
%% %%%%%%%%%%%%%%%%% Create a random test problem %%%%%%%%%%%%%%%%%%%%%%%%%
% Pick random matrix and vectors
m = 300;
n = 100;
A = 0.15 * rand(m, n);
b = randn(m, 1);

% Set the bounds
lb = zeros(n, 1) - 1;
ub = zeros(n, 1) + 1;

%% %%%%%%%%%%%%%%%% Solve the problem in different ways %%%%%%%%%%%%%%%%%%%
% Unconstraint
unconstraint_sol = (A' * A) \ (A' * b);
% Build-in MatLab function
quadprog_sol = quadprog(A' * A, -A' * b, [], [], [], [], lb, ub);
% ResQPASS (with quadprog)
[ResQPASS_sol, history] = ResQPASS_quadprog(A, b, lb, ub, "sol", quadprog_sol);
%% %%%%%%%%%%%%%%%% Determine related parameters to convergence %%%%%%%%%%%
% Calculate the error between the solutions
errorResQPASS = norm(ResQPASS_sol - unconstraint_sol);
errorQuadprog = norm(quadprog_sol - unconstraint_sol);

% Search the amount of active bounds which corresponds to the amount of 
% non-zero Lagrange multipliers.
tol = 1e-8;                     % Point from which to consider something non-zero

active_lower = zeros(n, 1);     % Active lower bounds
active_upper = zeros(n, 1);     % Active upper bounds

active_lower(history.lambda(:, end) > tol) = 1;
active_upper(history.mu(:, end) > tol) = 1; 

nr_bounds = sum(active_lower) + sum(active_upper);  % Count the amount of active bounds

% We see at iteration 32 that a value that exceeds the upperbound in the
% unconstraint solution thus not imply the value in the solution reaches
% the bound.
%% %%%%%%%%%%%%%%%%%%%%%%% Plot solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This figure aims to compare the solution of the ResQPASS-algorithm to
% that of the quadprog function and the unconstraint solution. The active
% bounds are also marked using Lagrange-multipliers.
figure("Name", "Comparisson of solutions", "NumberTitle", "off");
hold all;
markerSize=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all three solutions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ResQPASS_sol, "b");        % ResQPASS
plot(quadprog_sol, 'r--');      % quadprog
plot(unconstraint_sol, "g");    % Unconstrained

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the active bounds %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lowerbounds
idx = find(active_lower);
plot(idx, ResQPASS_sol(idx), "bo", ...
    "MarkerFaceColor", "b", "MarkerSize", markerSize);

% Upperbounds
idx = find(active_upper);
plot(idx, ResQPASS_sol(idx), "ro", ...
    "MarkerFaceColor", "r", "MarkerSize", markerSize);

%%%%%%%%%%%%%%%%%%%
% Plot the bounds %
%%%%%%%%%%%%%%%%%%%
line([0 n], [lb(1) lb(1)], "Color", "c");       % Lower bound
line([0 n], [ub(1) ub(1)], "Color", [245/250, 122/250, 126/250]);   % Upper bound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aditional figure settings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create axis labels
xlabel("$i$", "Interpreter", "latex");
ylabel("$x_i$", "Interpreter", "latex")

% Add legend
legend("ResQPASS", "quadprog", "Onbegrensd", ...
       "$\lambda_i \neq 0$", "$\mu_i \neq 0$", ...
       "Ondergrens", "Bovengrens", ...
       "Interpreter", "latex")
% Add title
title("Vergelijking tussen onbegrensde en begrensde oplossing")

% Save the figure
if saveFigure setFigParameters("onbegrensd-begrensd"); end

%% %%%%%%%%%%%%%%%%%%%% Plot residuals and objective %%%%%%%%%%%%%%%%%%%%%%
% This figure aims to show the convergence behaviour of the
% ResQPASS-algorithm with respect to the residuals and the objective
% function. 
figure("Name", "Residual and objective", "NumberTitle", "off");
%%%%%%%%%%%%%%%%%%%%%%%%%
% Examine the residuals %
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 2, 1);
semilogy(history.residuals);     % Plot the residuals
xline(nr_bounds, "--");          % Active bounds

% Extra plotting settings (labels, title and legend)
xlabel("Iteratie");             % x-label
ylabel("$\Vert A^T (Ax_k - b) - \lambda_k + \mu_k \Vert_2$", ...
    "Interpreter", "latex");    % y-label
title("Residu");                % Title
legend({"Residu", "Aantal actieve grenzen"}, "Location", "southwest")   % Legend

%%%%%%%%%%%%%%%%%%%%%%%%%
% Examine the objective %
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 2, 2);
semilogy(history.objective - norm(A*ResQPASS_sol - b)^2);   % Plot the objective
xline(nr_bounds, "--");                                     % Active bounds

% Extra plotting settings (title, labels, legend)
title("Objectief fout");                    % Title
xlabel("Iteratie");                         % x-label
ylabel("$\vert f(x^*) - f(x_k) \vert$", "Interpreter", "latex");                 % y-label
legend({"Objectief fout", "Aantal actieve grenzen"}, "Location", "southwest")    % Legend

% Save the results
if saveFigure setFigParameters("res_vs_obj"); end