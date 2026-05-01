% In this script, the way a preconditioner affects convergence of ResQPASS
% will be examined. The importance of a good first guess will also be
% pointed out.

% This corresponds to subsection 5.2.2 in the accompanying paper.

clear all;
close all;
clc;

% Determine whether the results should be saved
saveFigure = false;
%% %%%%%%%%%%%%%%%%%% Initialze the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the differentiation matrix A and pressure vector b
n = 50;
A = sparse(generate_2D(n));
b = 4*ones(n^2, 1);

% Make sure all methods work with the same parameters
x0 = zeros(n^2, 1);
tol = 1e-8;
maxiter = 2500;

% The boundaries for the solution
lb = zeros(n^2, 1);
ub = 0.1*ones(n^2, 1);

% Determine a preconditioner using an ilu-factorisation
setup.type = 'ilutp';
setup.droptol = 1e-3;

[L, U] = ilu(A' * A, setup);        % Perform ilu-factorisation
M = @(x) U\(L\x);                   % Preconditioner

% Choose diffrent initial guess
x1 = 1e-2*ones(n^2, 1);
%% %%%%%%%%%%%%%%%%%%%%%% Solve using LSQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lsqr_sol, lsqr_history] = lsqr_self(A, b, x0, maxiter, tol);

%% %%%%%%%%%%%%%%%%%%%%% Solve using ResQPass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notice that the problem is not solved right here. This is due to the fact
% that the problem is conditioned bad, meaning solving it takes a long
% time. So the solutions where saved and loaded later in the script.

% [resQPASS_sol, resQPASS_history] = ResQPASS_quadprog(A, b, lb, ub, ...
%     "maxIter", maxiter, ...
%     "x0", x0);
% save("PreConditioning_solutions.mat", "resQPASS_history", "resQPASS_sol", ...
%     "-append");

%% %%%%%%%%%%%%%%%%%%%%%%% Solve with better x0 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the same reason as before, the solution are saved and loaded later.

% [resQPASS_sol_init, resQPASS_history_init] = ResQPASS_altered(A, b, lb, ...
%     ub, "maxOuterIt", maxiter, "x0", x1, "maxOuterIt", 2000, ...
%     "displayResiduals", true);

% save("PreConditioning_solutions.mat", "resQPASS_history_init", ...
%     "resQPASS_sol_init", "-append");
%% %%%%%%%%%%%%%%%%%% Load the solutions to long problem %%%%%%%%%%%%%%%%%%
load("subsection_5_2_2_preConditioningSol.mat");

%% %%%%%%%%%%%%%%% Solve using LSQR with pre-conditioner %%%%%%%%%%%%%%%%%%
[lsqr_sol_precon, lsqr_history_precon] = lsqr_self(A, b, x0, maxiter, tol, M);

%% %%%%%%%%%%%% Solve using ResQPASS with pre-conditioner %%%%%%%%%%%%%%%%%
[resQPASS_sol_precon, resQPASS_history_precon] = ResQPASS_altered(A, b, lb, ub,...
    "maxOuterIt", maxiter, ...
    "maxInnerIt", 1e10, ...
    "M1", M, ...
    "x0", x0);

%% %%%%%%%%%%%%%%%%%%%%%%% Plot the solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name", "Solutions", "NumberTitle", "off");
x = linspace(0, 1, n+2);
y = linspace(0, 1, n+2);
[xgrid, ygrid] = meshgrid(x, y);

view_angle = [-37.5, 15];

%%%%%%%%%%%%%
% Unbounded %
%%%%%%%%%%%%%
subplot(1, 2, 1);
lsqr_nodes = add_boundary_vals(lsqr_sol);       % Convert solution to grid
surf(xgrid, ygrid, lsqr_nodes);                 % Plot the solution

% Axis labels
xlabel("x")
ylabel("y")
zlabel("u(x,y)")

view(view_angle);
zlim([0 0.25]);

title("Oplossing onbegrensd probleem")          % Title

%%%%%%%%%%%
% Bounded %
%%%%%%%%%%%
subplot(1, 2, 2);
resQPASS_nodes = add_boundary_vals(resQPASS_sol);   % Convert solution to grid
surf(xgrid, ygrid, resQPASS_nodes); hold all;       % Plot the solution

% Axis labels
xlabel("x")
ylabel("y")
zlabel("u(x,y)")

zlim([0 0.25]);
view(view_angle);

title("Oplossing begrensd probleem");               % Title

% Export image
if saveFigure, setFigParameters("Pressure-solution", 15.5, 6, [], 10); end

%% %%%%%%%%%%%%% Plot the solutions mid iteration %%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name", "Solutions at iteration 200", "NumberTitle", "off");

%%%%%%%%%%%%%%%%%%%%%%%%%%
% LSQR, unpreconditioned %
%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 2, 3);
nodes = add_boundary_vals(lsqr_history.solutions(:, 400));
surf(xgrid, ygrid, nodes);
xlabel("x");
ylabel("y");
zlabel("u(x, y)");
title({"Oplossing van LSQR", "op iteratie 400"});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ResQPASS unpreconditioned %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 2, 1);
nodes = add_boundary_vals(resQPASS_history.solutions(:, 1880));
surf(xgrid, ygrid, nodes);
xlabel("x");
ylabel("y");
zlabel("u(x, y)");
title({"Oplossing van ResQPASS", "op iteratie 1880"})

%%%%%%%%%%%%%%%%%%%%%%%%%
% LSQR - preconditioned %
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 2, 4);
nodes = add_boundary_vals(lsqr_history_precon.solutions(:, 4));
surf(xgrid, ygrid, nodes);
xlabel("x");
ylabel("y");
zlabel("u(x, y)");
title({"Oplossing van preconditioned",  "LSQR op iteratie 4"});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ResQPASS - preconditioned %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2, 2, 2);
nodes = add_boundary_vals(resQPASS_history_precon.solutions(:, 4));
surf(xgrid, ygrid, nodes);
xlabel("x");
ylabel("y");
zlabel("u(x, y)");
title({"Oplossing van preconditioned",  "ResQPASS op iteratie 4"});

% Export image
if saveFigure setFigParameters("Pressure-intermediat_sol", 15.5, 11, [], 10); end

%% %%%%%%%%%%%%%%%%%% Compare residuals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name", "Residuals", "NumberTitle", "off");

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without preconditioner %
%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = subplot(1, 2, 1);
semilogy(resQPASS_history.residuals); hold on    % Plot ResQPASS
semilogy(lsqr_history.residuals);               % Plot LSQR

xlabel("Iteratie");                             % x-label
ylabel("Residu");                               % y-label
title("Residuen zonder preconditioner")         % Title
legend("ResQPASS", "LSQR")                      % Legend

%%%%%%%%%%%%%%%%%%%%%%%
% With preconditioner %
%%%%%%%%%%%%%%%%%%%%%%%
ax2 = subplot(1, 2, 2);
semilogy(resQPASS_history_precon.residuals); hold on    % Plot ResQPASS
semilogy(lsqr_history_precon.residuals);                % Plot LSQR

xlabel("Iteratie");                                     % x-label
ylabel("Residu");                                       % y-label
title("Residuen met preconditioner")                    % Title
legend("ResQPASS", "LSQR")                              % Legend

% Export image
if saveFigure setFigParameters("Pressure-res", [], [], [], 10); end
%% %%%%%%%%%%%%% Growth of unpreconditioned solution - LSQR%%%%%%%%%%%%%%%%
figure("Name", "Growth unpreconditioned LSQR", "NumberTitle", "off")
steps = [10, 200, 400, 600, 1000, 1200]; % Iterations at which to display the results
view_angle = [-20, 50];

for k = 1:length(steps)
    subplot(3, 2, k);
    nodes = add_boundary_vals(lsqr_history.solutions(:, steps(k)));
    % Plot the growth of the lsqr solution at specified iterations
    surf(xgrid, ygrid, nodes);
    title(sprintf("Oplossing van LSQR op iteratie %d", steps(k)));
    view(view_angle)
    xlabel("x")
    ylabel("y")
    zlabel("u(x, y)")
end

% Export image
if saveFigure setFigParameters("Pressure-LSQR_growth", 15.5, 20, [], 10); end

%% %%%%%%%%%% Growth of unpreconditioned solution - ResQPASS %%%%%%%%%%%%%%
steps = [400, 1000, 1600, 1800, 1900, 2400]; % Iterations at which to display the results
figure("Name", "Growth unpreconditioned ResQPASS", "NumberTitle", "off");

for k = 1:length(steps)
    subplot(3, 2, k);
    nodes = add_boundary_vals(resQPASS_history.solutions(:, steps(k)));
    % Plot the growth of the preconditioned solution at specified steps
    surf(xgrid, ygrid, nodes);
    title({"Oplossing van ResQPASS", sprintf("op iteratie %d", steps(k))});
    view(view_angle)
    xlabel("x")
    ylabel("y")
    zlabel("u(x, y)")
end

% Export image
if saveFigure setFigParameters("Pressure-ResQPASS_growth", 15.5, 20, [], 10); end

%% %%%%%%%%%%%% Growth of preconditioned solution - ResQPASS %%%%%%%%%%%%%%
steps = round(linspace(1, 15, 6));

figure("Name", "Growth predonditioned ResQPASS", "NumberTitle", "off");

for k = 1:length(steps)
    subplot(3, 2, k);
    nodes = add_boundary_vals(resQPASS_history_precon.solutions(:, steps(k)));
    % Plot the growth of the preconditioned solution at specified steps
    surf(xgrid, ygrid, nodes);
    title({"Oplossing van preconditioned",  sprintf("ResQPASS op iteratie %d", steps(k))});
    view(view_angle)
    xlabel("x")
    ylabel("y")
    zlabel("u(x, y)")
end

% Export image
if saveFigure setFigParameters("Pressure-ResQPASS_growth_precondition", 15.5, 20, [], 10); end

%% %%%%%%%%%%%%% Growth of preconditioned solution - LSQR  %%%%%%%%%%%%%%%%
steps = round(linspace(1, 12, 6));
figure("Name", "Growth preconditioned LSQR", "NumberTitle", "off");

for k = 1:length(steps)
    subplot(3, 2, k);
    nodes = add_boundary_vals(lsqr_history_precon.solutions(:, steps(k)));
    % Plot the growth of the preconditioned solution at specified steps
    surf(xgrid, ygrid, nodes);
    title({"Oplossing van preconditioned",  sprintf("LSQR op iteratie %d", steps(k))});
    view(view_angle)
    xlabel("x")
    ylabel("y")
    zlabel("u(x, y)")
end

% Export image
if saveFigure setFigParameters("Pressure-LSQR_growth_preconditioning", 15.5, 20, [], 10); end

%% %%%%%%%%%%%%%%%%%%%%%%%%% Basis vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this figure, a heatmap will be shown corresponding to the basis
% vectors of V.
figure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without preconditioning %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 2, 1)
V_trim = resQPASS_history.V(:, 1:2460); % Shorten V to cut it into blocks
% Dimension of the blocks
blockSizeRow = 4;
blockSizeCol = 4;

% Number of blocks
nRows = size(V_trim,1)/blockSizeRow;
nCols = size(V_trim,2)/blockSizeCol;

% Reshape V to get the blocks
V_reshaped = reshape(V_trim, blockSizeRow, nRows, blockSizeCol, nCols);

% Take the average of all blocks to obtain lower resolution
V_small = squeeze(mean(mean(V_reshaped,1),3));

% Flatten the outliers to make the diffrences more visable
z = V_small(:);
low = prctile(z, 1);
high = prctile(z, 99);

% Plot the results
imagesc([1, 2460], [1, 2500], V_small);
set(gca,'YDir','normal')                    % Make ticks go low to high
caxis([low high]);                          % Add the right indecis
xlabel("Iteratie j");                 % x-label
ylabel("Residu component (r_j)_i");   % y-label

title({"Heatmap van basisvectoren", "zonder preconditioning (V_{ij})"}); % Title

ax = gca;
ti = ax.TightInset;
ax.Position = [ ...
    ax.Position(1) + ti(1), ...
    ax.Position(2) + ti(2), ...
    ax.Position(3) - ti(1) - ti(3), ...
    ax.Position(4) - ti(2) - ti(4) ...
];
%%%%%%%%%%%%%%%%%%%%%%%%
% With preconditioning %
%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 2, 2)

% Flatten for more clear diffrences
z = resQPASS_history_precon.V(:);
imagesc(resQPASS_history_precon.V);         % Heatmap of the basis vectors
set(gca,'YDir','normal')                    % Make ticks go low to high
xlabel("Iteratie j");                 % x-label
ylabel("Residu component (r_j)_i");   % y-label

title({"Heatmap van basisvectoren", "met preconditioning (V_{ij})"})    % Title

ax = gca;
ti = ax.TightInset;
ax.Position = [ ...
    ax.Position(1) + ti(1), ...
    ax.Position(2) + ti(2), ...
    ax.Position(3) - ti(1) - ti(3), ...
    ax.Position(4) - ti(2) - ti(4) ...
];

if saveFigure, exportgraphics(fig, "Images/Pressure-Heatmap.pdf","ContentType", ...
"image", "Resolution", 300); end
%% %%%%%%%%%%%%%%%%%%%%%% Explain assymetry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,D] = eigs(inv(L*U) * (A'*A), 5);
% V = V(:);
imagesc(real(V'))
yticks(1:5);
yticklabels({
    "Eigenvector 1", ...
    "Eigenvector 2", ...
    "Eigenvector 3", ...
    "Eigenvector 4", ... 
    "Eigenvector 5"
});
xlabel("Component index");
set(gca,'YDir','normal');

title("5 Grootste eigenvectoren van (LU)^{-1}(A^T A)");
if saveFigure, setFigParameters("Pressure-assymetry", 15.5, 3, [], 10); end

%% %%%%%%%%%%%%%%%%%%%%%%%% Active bounds  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name", "Active bounds")
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without preconditioner %
%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)

mu_with_bounds = add_boundary_vals(resQPASS_history.mu(:,2400));
mask  = mu_with_bounds > 1e-8;

imagesc(x, y, resQPASS_nodes);
set(gca,'YDir','normal');

axis equal
hold on

[row, col] = find(mask);        % places where mask is true

scatter(x(row), y(col), 10, 'r', 'filled');

xlim([0 1]);
ylim([0 1]);
legend("Actieve grenzen");
xlabel("x");
ylabel("y");
title({"Actieve grenzen", "zonder preconditionering"});

%%%%%%%%%%%%%%%%%%%%%%%%%%
% With preconditioner %
%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)

mu_with_bounds = add_boundary_vals(resQPASS_history_precon.mu);
mask  = mu_with_bounds > 1e-8;

resQPASS_precon_nodes = add_boundary_vals(resQPASS_sol_precon);

imagesc(x, y, resQPASS_precon_nodes);
set(gca,'YDir','normal');

axis equal
hold on

[row, col] = find(mask);      % places where mask is true

scatter(x(row), y(col), 10, 'r', 'filled');

xlim([0 1]);
ylim([0 1]);
legend("Actieve grenzen");
xlabel("x");
ylabel("y");
title({"Actieve grenzen", "met preconditionering"});

if saveFigure, setFigParameters("Pressure-active_bounds", [], [], [], 15); end
%% %%%%%%%%%%%%%%%%% Extra hypothesis: initial guess %%%%%%%%%%%%%%%%%%%%%%
% The first thing that was noticed is that the choice of x0 can greatly
% affect the convergence. When chosing x0 = 1e-3 instead of x0 = 0, a close
% solution is found much sooner.
figure;
%%%%%%%%%%
% x0 = 0 %
%%%%%%%%%%
subplot(1, 2, 1)
nodes = add_boundary_vals(resQPASS_history.solutions(:, 1000));
surf(xgrid, ygrid, nodes);
title({"Oplossing op iteratie 1000", "met x0 = 0"})
%%%%%%%%%%%%%
% x0 = 1e-3 %
%%%%%%%%%%%%%
subplot(1, 2, 2)
nodes = add_boundary_vals(resQPASS_history_init.solutions(:, 1000));
surf(xgrid, ygrid, nodes);
title({"Oplossing op iteratie 1000", "met x0 = 0.001"})
if saveFigure, setFigParameters("Pressure-diffrent_initVal", 15.5, 6); end

%% %%%%%%%%%%%%%% Effect halving characteristic length %%%%%%%%%%%%%%%%%%%%
n = 25;
A = generate_2D(n);
b = 4*ones(n^2, 1);

% Make sure all methods work with the same parameters
x0 = zeros(n^2, 1);
tol = 1e-8;
maxiter = 2500;

% The boundaries for the solution
lb = zeros(n^2, 1);
ub = 0.1*ones(n^2, 1);

% Determine a preconditioner using an ilu-factorisation
setup.type = 'ilutp';
setup.droptol = 1e-3;

[L, U] = ilu(A' * A, setup);        % Perform ilu-factorisation
M = @(x) U\(L\x);                   % Preconditioner

[sol_n25, history_n25] = ResQPASS_altered(A, b, lb, ub, ...
    "maxOuterIt", maxiter, ...
    "maxInnerIt", 200, ...
    "M1", M, ...
    "x0", x0, ...
    "displayRes", true);

n = 101;
A = generate_2D(n);
b = 4*ones(n^2, 1);

% Make sure all methods work with the same parameters
x0 = zeros(n^2, 1);
tol = 1e-8;
maxiter = 2500;

% The boundaries for the solution
lb = zeros(n^2, 1);
ub = 0.1*ones(n^2, 1);

% Determine a preconditioner using an ilu-factorisation
setup.type = 'ilutp';
setup.droptol = 1e-3;

[L, U] = ilu(A' * A, setup);        % Perform ilu-factorisation
M = @(x) U\(L\x);                   % Preconditioner

[sol_n101, history_n101] = ResQPASS_altered(A, b, lb, ub, ...
    "maxOuterIt", maxiter, ...
    "maxInnerIt", 200, ...
    "M1", M, ...
    "x0", x0, ...
    "displayRes", true);
%%
% Extract Lagrange multipliers for different grid sizes
mu25 = history_n25.mu;
lambda25 = history_n25.lambda;

% Count number of active constraints (thresholded)
active_bounds_n25 = sum(mu25 > 1e-8) + sum(lambda25 > 1e-8);

mu50 = resQPASS_history_precon.mu;
lambda50 = resQPASS_history_precon.lambda;
active_bounds_n50 = sum(mu50 > 1e-8) + sum(lambda50 > 1e-8);

mu101 = history_n101.mu;
lambda101 = history_n101.lambda;
active_bounds_n101 = sum(mu101 > 1e-8) + sum(lambda101 > 1e-8);

% Display number of active constraints for each discretization
disp("Amount of active bounds n = 25: " + active_bounds_n25);
disp("Amount of active bounds n = 50: " + active_bounds_n50);
disp("Amount of active bounds n = 101: " + active_bounds_n101);
disp("Amount of iterations n = 25: " + history_n25.outerIt);
disp("Amount of iterations n = 50: " + resQPASS_history_precon.outerIt);
disp("Amount of iterations n = 101: " + history_n101.outerIt);
figure;

%%%%%%%%%%
% n = 25 %
%%%%%%%%%%
subplot(1,3,1); hold on;

n = 25;

% Generate grid coordinates (including boundary points)
x = linspace(0,1,n+2);
y = linspace(0,1,n+2);

% Reconstruct solution including boundary values
nodes = add_boundary_vals(sol_n25);

% Plot solution as background image
imagesc(x, y, nodes);

% Find indices of active constraints (mu > tolerance)
[row,col] = find(add_boundary_vals(mu25) > 1e-8);

% Overlay active constraint locations
scatter(x(col), y(row), 2, 'r', 'filled');

% Formatting
axis equal;
xlim([0 1]); ylim([0 1]);
xlabel("x")
ylabel("y")
title({"Actieve grenzen", "met n = 25"});


%%%%%%%%%%
% n = 50 %
%%%%%%%%%%
subplot(1,3,2); hold on;

n = 50;

% Generate grid coordinates (including boundary points)
x = linspace(0,1,n+2);
y = linspace(0,1,n+2);

% Reconstruct solution including boundary values
nodes = add_boundary_vals(resQPASS_sol_precon);

% Plot solution as background image
imagesc(x, y, nodes);

% Find indices of active constraints (mu > tolerance)
[row,col] = find(add_boundary_vals(mu50) > 1e-8);

% Overlay active constraint locations
scatter(x(col), y(row), 2, 'b', 'filled');

% Formatting
axis equal;
xlim([0 1]); ylim([0 1]);
xlabel("x")
ylabel("y")
title({"Actieve grenzen", "met n = 50"});


%%%%%%%%%%%
% n = 101 %
%%%%%%%%%%%
subplot(1,3,3); hold on;

n = 101;

% Generate grid coordinates (including boundary points)
x = linspace(0,1,n+2);
y = linspace(0,1,n+2);

% Reconstruct solution including boundary values
nodes = add_boundary_vals(sol_n101);

% Plot solution as background image
imagesc(x, y, nodes);

% Find indices of active constraints (mu > tolerance)
[row,col] = find(add_boundary_vals(mu101) > 1e-8);

% Overlay active constraint locations
scatter(x(col), y(row), 2, 'g', 'filled');

% Formatting
axis equal;
xlim([0 1]); ylim([0 1]);
xlabel("x")
ylabel("y")
title({"Actieve grenzen", "met n = 101"});

% Save image
if saveFigure, setFigParameters("Pressure-active_bounds_n", [], [], [], 10); end

%% Functions to generated the Laplacian matrix and add boundary points.
% Generate Laplacian approximation
function A = generate_1D(n)
    A = (n + 1)^2 * spdiags([-ones(n,1)  2*ones(n,1)  -ones(n,1) ], [-1 0 1], n, n);
end

function A = generate_2D(n)
    A1D = generate_1D(n);
    A = kron(speye(n), A1D) + kron(A1D, speye(n));
end

% Add booundary points
function nodes = add_boundary_vals(vec)
    n = sqrt(size(vec, 1));
    nodes = zeros(n+2, n+2);
    nodes(2:end-1, 2:end-1) = reshape(vec, n, n);
end
