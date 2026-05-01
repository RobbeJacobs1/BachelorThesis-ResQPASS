% This script researches the effect of the amount of active bounds to the
% convergence behaviour of the ResQPASS-algorithm.

% This corresponds to Example 4.2.1 (Figure 4.3, p. 24) in the accompanying paper.
clear all; close all; clc;
rng(42);

% Determine whether the results should be saved.
saveFigure = true;
%% %%%%%%%%%%%%%%%%%%%%% Setup up the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 600;
n = 400;

% Generate a matrix A with a good condition numer
U = orth(randn(m, n));
V = orth(randn(n,n));
kappa = 2;
S = linspace(1, 1/kappa, n);
A = U * diag(S) * V;

lb = -Inf * ones(n, 1);     % No lower bound
ub = ones(n, 1);            % Upper bound at 1

x = -0.5 + 1 * rand(n, 1);  % Random unconstraint solution well between the bounds

%% %%%%%%%%Solve the problem for diffrent number of bounds %%%%%%%%%%%%%%%%
bounds = 0:10:100;                  % Determine number of active bounds
offsets = 1.5 + 0.5*rand(n, 1);     % Randomly select the values of the solutions at the active bounds
idx = randperm(n, 400)';            % Select which bounds should become active

residuals = zeros(400, size(bounds, 2));
for k = 1:size(bounds, 2)
    nr_bounds = bounds(k);                              
    x(idx(1:nr_bounds)) = offsets(idx(1:nr_bounds));    % Put certain parts of the solution behind the bouds
    b = A * x;                                          % Make x the actual solution
    [sol, history] = ResQPASS_quadprog(A, b, lb, ub, "sol", x); % Solve the problem
    residuals(1:history.iters, k) = history.residuals;   % Save the residuals
end

%% %%%%%%%%%%%%%%%%%%%%%% Plot the residuals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on
yscale log

nr_tests = size(residuals,2);               % Amount of tests
residuals_plots = gobjects(nr_tests,1);     % to store semilogy handles
for k = 1:nr_tests
    residuals_plots(k) = semilogy(residuals(:,k), ...   % Plot the residuals
        'DisplayName', sprintf('%d grenzen', bounds(k)));
end

% Add axis labels
xlabel('Iteratie');
ylabel('Residu');

% Add vertical lines at "nr_bounds" iteration
for k = 1:nr_tests
    % Pick the color of the corresponding curve
    c = residuals_plots(k).Color;
    
    % Draw vertical line at iteration = bounds(k)
    % (or any other x-position representing bound identification)
    xline(bounds(k), '--', sprintf('%d', bounds(k)), ...
          'Color', c, 'LabelOrientation','horizontal', ...
          'LabelVerticalAlignment','bottom', ...
          "HandleVisibility", "off");
end

% Add legend and title
legend('show');
title("Vergelijking tussen het residu en het aantal actieve grenzen");

% Save results
if saveFigure setFigParameters("Bounds"); end
