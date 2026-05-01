% In this script, the effect of halving the characteristic length of a mesh
% has when the mesh is squeezed using box constraints. The rectangle
% example will not be covered further as it gave no interesting results
% diffrent from the other two cases.

% This corresponds to section 6.2 in the accompanying paper.

close all; clear; clc

% Determine whether the figures should be saved
saveFigure = false;
%% %%%%%%%%%%%%%%%%%%%%%% Show solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
% Load results %
%%%%%%%%%%%%%%%%

% Circle
circle = load("Create_2D_meshes/Results/circle_0.05.mat");
circle.name = "schijf";

% H
H = load("Create_2D_meshes/Results/H_0.05.mat");
H.name = "H";

% Rectangle
rect = load("Create_2D_meshes/Results/rectangle_0.05.mat");
rect.name = "rechthoek";

%%%%%%%%%%%%%%%%
% Plot results %
%%%%%%%%%%%%%%%%
% Circle
figure("Name", "Circle", "NumberTitle", "off");
ax1 = axes();
show_solution(circle, ax1);
if saveFigure, setFigParameters("deformation-circle_sol"); end

% H
figure("Name", "H", "NumberTitle", "off");
ax2 = axes();
show_solution(H, ax2);
if saveFigure, setFigParameters("deformation-H_sol"); end

% Rectangle
figure("Name", "Rectangle", "NumberTitle", "off");
ax3 = axes();
show_solution(rect, ax3);
if saveFigure, setFigParameters("deformation-rectangle_sol"); end
%% %%%%%%%%%%%%%%%%%%% Compare residual to refinement %%%%%%%%%%%%%%%%%%%%%
h = [0.05, 0.07, 0.1, 0.14, 0.2]; % Select refinements to compare

figure("Name", "Residuals of circle", "NumberTitle", "off")
hold on;
yscale("log");

for k = 1:length(h)
    msh = load("Create_2D_meshes/Results/circle_" + ...
        sprintf("%.2f", h(k)) + ".mat");                            % Load mesh
    semilogy(msh.history.residuals, "DisplayName", "h = " + h(k));  % Plot the results
end
title("Residuen schijf")                                            % Title
legend()                                                            % Legend
if saveFigure, setFigParameters("deformation-res_circle"); end      % Save results

figure("Name", "Residuals of H", "NumberTitle", "off")
hold on;
yscale("log");
for k = 1:length(h)
    msh = load("Create_2D_meshes/Results/H_" + ...
        sprintf("%.2f", h(k)) + ".mat");                            % Load mesh
    semilogy(msh.history.residuals, "DisplayName", "h = " + h(k));  % Plot the results
end
title("Residuen H")                                                 % Title
legend()                                                            % Legend
if saveFigure, setFigParameters("deformation-res_H"); end           % Save results
%% %%%%%%%%%%%%%%%%%%%% Compare number of active bounds %%%%%%%%%%%%%%%%%%%
figure("Name", "Amount of active bounds, circle", "NumberTitle", "off");

%%%%%%%%%%
% Circle %
%%%%%%%%%%
subplot(1, 2, 1);
hold on;
h = 0.05:0.01:0.20;                                         % Refinements
active_bounds = zeros(length(h), 1);
for k = 1:length(h)
    msh = load("Create_2D_meshes/Results/circle_" + ...
        sprintf("%.2f", h(k)) + ".mat");                    % Load mesh
    tol = 1e-8;                                             % What to consider zero
    active_lower = sum(msh.history.lambda(:, end) > tol);   % Amount of active lower bounds
    active_upper = sum(msh.history.mu(:, end) > tol);       % Amount of active upper bounds
    active_bounds(k) = active_lower + active_upper;         % Amount of active bounds
end
scatter(h, active_bounds);                                  % Plot the results

% Plotting settings
xlabel("Karakteristieke lengte")                            % x-label
ylabel("Aantal actieve grenzen")                            % y-label
title({"Effect karakteristieke lengte, ", ...
    "aantal actieve grenzen (schijf)"})                     % Title

%%%%%
% H %
%%%%%
subplot(1, 2, 2)
hold on;
h = 0.05:0.01:0.20;                                         % Refinements
active_bounds = zeros(length(h), 1);            
for k = 1:length(h)
    msh = load("Create_2D_meshes/Results/H_" + ...
        sprintf("%.2f", h(k)) + ".mat");                    % Load mesh
    tol = 1e-8;                                             % Tolerance to determine zero
    active_lower = sum(msh.history.lambda(:, end) > tol);   % Active lower bounds
    active_upper = sum(msh.history.mu(:, end) > tol);       % Active upper bounds
    active_bounds(k) = active_lower + active_upper;         % Amount of active bounds
end
scatter(h, active_bounds);                                  % Plot the results

% Plotting settings
xlabel("Karakteristieke lengte")                            % x-label
ylabel("Aantal actieve grenzen")                            % y-label
title({"Effect karakteristieke lengte,",  ...
    "aantal actieve grenzen (H)"})                          % Title
if saveFigure, setFigParameters("deformation-nrBounds"); end  % Save results

%% %%%%%%%%%%%%%%%%% Show effect of doubling nodes %%%%%%%%%%%%%%%%%%%%%%%%
% This figure shows zoomed in versions of the circle to examine the effect
% of halving the characteristic length.

figure("Name", "Effect of doubling h (edge)", "NumberTitle","off");

% Load the meshes together with their properties
circle5 = load("Create_2D_meshes/Results/circle_0.05.mat");
triangles5 = circle5.triangles;
sol5 = circle5.sol_nodes;

circle10 = load("Create_2D_meshes/Results/circle_0.10.mat");
sol10 = circle10.sol_nodes;
triangles10 = circle10.triangles;
nodes10 = circle10.nodes;

% ------------------ Active bounds --------------
lambdax = circle5.history.lambda(1:2:end, end);    % Lagrange multipliers
lambday = circle5.history.lambda(2:2:end, end);

mux = circle5.history.mu(1:2:end, end);
muy = circle5.history.mu(2:2:end, end);

% Determine the active bounds
activeLowerX = lambdax > tol;
activeUpperX = mux > tol;
activeLowerY = lambday > tol;
activeUpperY = muy > tol;

% ------------------ Zoom area ------------------
zoom_x = [-0.9, -0.5];
zoom_y = [0, 0.4];

% ------------------ Zoomed plot ------------------
ax1 = axes('Position',[0.1 0.1 0.75 0.8]);
hold on;

% Plot the zoomed plot
triplot(triangles5, sol5(:, 1), sol5(:, 2), "r");
triplot(triangles10, sol10(:, 1), sol10(:, 2), "b", "LineWidth", 1.2);

% Add the boundary points
scatter(sol5(activeLowerX,1), sol5(activeLowerX,2), "r");
scatter(sol5(activeUpperX,1), sol5(activeUpperX,2), "r");
scatter(sol5(activeLowerY,1), sol5(activeLowerY,2), "r");
scatter(sol5(activeUpperY,1), sol5(activeUpperY,2), "r");

axis equal;

% apply the zoom
xlim(zoom_x);
ylim(zoom_y);

% Title
title({"Effect h verdubbelen", "op actieve grenzen"});

% ------------------ Full plot ------------------
ax2 = axes('Position',[0.45 0.58 0.28 0.28]);
hold on;

% Plot the whole result
triplot(triangles5, sol5(:, 1), sol5(:, 2), "r");
triplot(triangles10, sol10(:, 1), sol10(:, 2), "b");
% Remove the ticks
set(gca, 'XTick', [], 'YTick', []);
axis equal;
% Make the limits square
xlim([-1 1]);
ylim([-1 1]);

% -------------------- Zoom box ---------------------
% Add a box around the zoomed area
axes(ax2);
rectangle('Position', ...
    [zoom_x(1), zoom_y(1), diff(zoom_x), diff(zoom_y)], ...
    'EdgeColor','k','LineWidth',1.5);

% ------------------ ax2 on top ------------------
uistack(ax2, 'top');

% Save the results
if saveFigure, setFigParameters("deformation-effect_doubling_h"); end

%% %%%%%%%%% The solutions tries to preserve original distances %%%%%%%%%%%
figure("Name", "Effect of doubling h (internal)", "NumberTitle","off");

% ------------------ Zoom area -----------------
zoom_x = [0, 0.4];
zoom_y = [0, 0.4];

% ------------------ Zoomed plot ---------------
ax1 = axes('Position',[0.1 0.1 0.75 0.8]);
hold on;

% Plot the meshes
triplot(triangles10, nodes10(:, 1), nodes10(:, 2), "g");
triplot(triangles10, sol10(:, 1), sol10(:, 2), "r");

axis equal;

% Apply the zoom
xlim(zoom_x);
ylim(zoom_y);

% Title
title("Effect vervorming interne structuur");

% ------------------ Full plot ------------------
ax2 = axes('Position',[0.45 0.58 0.28 0.28]);
hold on;

% Plot the results
triplot(triangles10, nodes10(:, 1), nodes10(:, 2), "g");
triplot(triangles10, sol10(:, 1), sol10(:, 2), "r");

% Remove ticks
set(gca, 'XTick', [], 'YTick', []);
axis equal;
% Make the limits square
xlim([-1 1]);
ylim([-1 1]);


% -------------------- Zoom box ---------------------
% Add a box around the zoomed area
rectangle('Position', ...
    [zoom_x(1), zoom_y(1), diff(zoom_x), diff(zoom_y)], ...
    'EdgeColor','k','LineWidth',1.5);

% ------------------ ax2 on top --------------------
uistack(ax2, 'top');
if saveFigure, setFigParameters("deformation-internal_structure"); end

%% %%%%%%%%%%%%%%%%%%% Compare with preconditioning %%%%%%%%%%%%%%%%%%%%%%%
% Load the solutions without preconditioning
circle_5 = load("Create_2D_meshes/Results/circle_0.05.mat");
circle_10 = load("Create_2D_meshes/Results/circle_0.10.mat");
circle_20 = load("Create_2D_meshes/Results/circle_0.20.mat");
circles = {circle_5, circle_10, circle_20};             % Container for all results

% Load the solutions with preconditioning
pcircle_5 = load("Create_2D_meshes/Results/circle_0.05_precond.mat");
pcircle_10 = load("Create_2D_meshes/Results/circle_0.10_precond.mat");
pcircle_20 = load("Create_2D_meshes/Results/circle_0.20_precond.mat");
pcircles = {pcircle_5, pcircle_10, pcircle_20};         % Container for all results

figure("Name", "Residuals with preconditioning", "NumberTitle", "off");
hold on;
yscale("log");
c = get(gca, 'ColorOrder');                     % Get color palet

alpha = 0.7;                                    % Amount of darkening

for k = 1:length(circles)
    c_dark = alpha*c(k, :);                     % Darker color variant

    % Without preconditioning
    semilogy(circles{k}.history.residuals, "Color", c(k,:));
    
    % With preconditioning (darker)
    semilogy(pcircles{k}.history.residuals, "Color", c_dark);
end

% Legend
legend("h = 0.05", "h = 0.10", "h = 0.20");

xlabel("Iteratie");                             % x-label
ylabel("Residu");                               % y-label
title("Effect preconditioner")                  % Title
if saveFigure, setFigParameters("deformation-precon"); end  % Save result


function ax = show_solution(mesh, ax)
% This function gives a way to present a mesh with its boundary points.
% ----------------- Get information about the mesh ------------------------
    nodes = mesh.nodes;                             % nodes
    solution = mesh.sol_nodes;                      % Solution
    triangles = mesh.triangles;                     % Triangles

    lambdax = mesh.history.lambda(1:2:end, end);    % Lagrange multipliers
    lambday = mesh.history.lambda(2:2:end, end);
    
    mux = mesh.history.mu(1:2:end, end);
    muy = mesh.history.mu(2:2:end, end);

    lb_x = mesh.lb_x;                               % Boundaries
    ub_x = mesh.ub_x;
    lb_y = mesh.lb_y;
    ub_y = mesh.ub_y;


    axes(ax);                                      % Select axis as current
    hold(ax, "on")

% ----------------- Original and deformed mesh ----------------------------
    triplot(triangles, nodes(:, 1), ...
        nodes(:, 2), "g-");                         % original mesh
    triplot(triangles, solution(:, 1), ...
        solution(:, 2), "r-")                       % Deformed mesh
    

% ----------------- Determine active boundaries ---------------------------
    tol = 1e-8;
    activeLowerX = lambdax > tol;
    activeUpperX = mux > tol;
    activeLowerY = lambday > tol;
    activeUpperY = muy > tol;
    
% ------------------------- Boundary points--------------------------------
    scatter(solution(activeLowerX, 1), ...          % Active left bounds
        solution(activeLowerX, 2), ...
        "b");
    
    scatter(solution(activeUpperX, 1), ...          % Active right bounds
        solution(activeUpperX, 2), ...
        "b");
    
    scatter(solution(activeLowerY, 1), ...          % Active lower bounds
        solution(activeLowerY, 2), ...
        "b");
    
    scatter(solution(activeUpperY, 1), ...          % Active upper bounds
        solution(activeUpperY, 2), ...
        "b");

% -------------------------- Boundary lines -------------------------------
    x = [lb_x, ub_x, ub_x, lb_x, lb_x];
    y = [lb_y, lb_y, ub_y, ub_y, lb_y];
    line(x, y, "Color", [0 0 0 0.3], "LineWidth", 2);

   % -------------------------- Set plot bounds ---------------------------
    minx = min([nodes(:, 1); solution(:, 1)]);     % Lowest x-value
    maxx = max([nodes(:, 1); solution(:, 1)]);     % Highest x-value
    miny = min([nodes(:, 2); solution(:, 2)]);     % Lowest y-value
    maxy = max([nodes(:, 2); solution(:, 2)]);     % Highest y-value
    
    % Change the plotting area
    xlim([minx - 0.25, maxx + 0.25]);
    ylim([miny - 0.25, maxy + 0.25]);

% ----------------------- Additional plot settings ------------------------
    axis square;
    legend("Originele mesh", "Vervormde mesh", "Actieve grenzen");
    title("Oplossing voor " + mesh.name);
end