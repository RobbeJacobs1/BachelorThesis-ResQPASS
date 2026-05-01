% In this script, the effect of the condition number on the convergence of
% the ResQPASS-algorithm will be researched and compared to that of LSQR.

% This corresponds to Example 4.2.2 (Figure 4.4, p. 25) in the accompanying paper.
clear all; close all; clc;
rng(42);

% Determine whether the results should be saved.
saveFigure = false;
%% %%%%%%%%%%%%%%%%%%%%% Setup up the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 600;
n = 400;

% Generate a matrix A with a set condition number using eigenvalue
% decompositions.
U = orth(randn(m, n));
V = orth(randn(n,n));

lb = -ones(n, 1);       % No lower bound
ub = ones(n, 1);        % Upper bound at 1
x = randn(n, 1);        % Choose the solution
x0 = zeros(n, 1);       % Set the initial guess for better comparisson.   

%% %%%%%% Solve the problem for diffrent Condition numbers %%%%%%%%%%%%%%%%
residuals_ResQPASS = zeros(400, 6);
residuals_lsqr = zeros(400, 6);
for k = 1:6
    kappa = 2^(k);                      % Set the condition number
    S = logspace(0, -log10(kappa), n);  % Choose eigenvalues logarithmically distributed
    A = U * diag(S) * V;                % Create A with condition number kappa
    b = A * x;                          % Fix be to make x the solution

    % Solve using ResQPASS
    [sol, ResQPASS_history] = ResQPASS_quadprog(A, b, lb, ub, "x0", x0);
    residuals_ResQPASS(1:ResQPASS_history.iters, k) = ...
        ResQPASS_history.residuals/ResQPASS_history.residuals(1);  % Extract residuals from history

    % Solve using LSQR
    [~, lsqr_history] = lsqr_self(A, b, x0, 1000, 1e-8);
    residuals_lsqr(1: lsqr_history.iters, k) = ...                 % Get residuals
        lsqr_history.residuals/lsqr_history.residuals(1);
end
%% %%%%%%%%%%%%%%%%%%%%%%% Plot the residuals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
cond_numbers = 2.^(1:6);

%%%%%%%%%%%%
% ResQPASS %
%%%%%%%%%%%%
ax1 = subplot(1,2,1);
hold on
yscale log

nr_tests = size(residuals_ResQPASS,2);
for k = 1:nr_tests
    semilogy(residuals_ResQPASS(:,k), ...       % Plot the residuals for each condition number
        'DisplayName', sprintf('$\\kappa(A) = %d$', cond_numbers(k)));
end

% Add a vertical line for the number of active bounds
xline(nnz(x > 1) + nnz(x < -1), "k--", "DisplayName", ...
    "$\vert \mathcal{W} \vert$", "Interpreter", "latex")

% Set x- and y-labels
xlabel('Iteratie $k$', "Interpreter", "latex");
ylabel("$\displaystyle \frac{\Vert r_k \Vert}{\Vert r_0 \Vert}$", "Interpreter", "latex", ...
    "Rotation", 0);

% Add legend and title
legend('show', 'Interpreter', 'latex');
title("Residuen van ResQPASS");

%%%%%%%%
% LSQR %
%%%%%%%%
ax2 = subplot(1,2,2);
hold on
yscale log

nr_tests = size(residuals_lsqr,2);
for k = 1:nr_tests
    semilogy(residuals_lsqr(:,k), ...           % Plot the residuals
        'DisplayName', sprintf('$\\kappa(A) = %d$', cond_numbers(k)));
end

% Add labels, legend and title
xlabel('Iteratie $k$', "Interpreter", "latex");
ylabel("$\frac{\Vert r_k \Vert}{\Vert r_0 \Vert}$", "Interpreter", "latex", ...
    "Rotation", 0);
legend('show', 'Interpreter', 'latex');
title("Residuen van LSQR");

% Force identical limits
linkaxes([ax1 ax2],'xy');
ylim([1e-8, Inf])
xlim([-Inf, 400])

% Save the figure
if saveFigure setFigParameters("conditiegetal"); end