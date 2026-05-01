% This script shows the effect of halving the characteristic length of a
% mesh for contact problems.

% This corresponds to section 6.3 in the accompanying paper.
close all; clear; clc

% Determine whether the figures should be saved
saveFigure = false;
%% %%%%%%%%%%%%%%%%%%%%%% Load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sphere = load("Create_3D_meshes/Results/sphere_0.10.mat");

%%%%%%%%%%%%%%%%%%%%% Show the test results %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name", "Result", "NumberTitle", "off");
hold on;
view(3);
% ----------------- Get information about the mesh ------------------------
nodes = sphere.nodes;                             % nodes
solution = sphere.sol_nodes;                      % Solution
triangles = sphere.outer_faces;                   % Triangles

lambdax = sphere.history.lambda(1:3:end, end);    % Lagrange multipliers
lambday = sphere.history.lambda(2:3:end, end);
lambdaz = sphere.history.lambda(3:3:end, end);

mux = sphere.history.mu(1:3:end, end);
muy = sphere.history.mu(2:3:end, end);
muz = sphere.history.mu(3:3:end, end);

% ----------------- Original and deformed mesh ----------------------------
trisurf(triangles, nodes(:, 1), ...              % original mesh
    nodes(:, 2), nodes(:, 3), ...
    "FaceColor", "none", "EdgeColor", "red"); 
trisurf(triangles, solution(:, 1), ...
    solution(:, 2), solution(:, 3))              % Deformed mesh
colormap(summer)

% ----------------- Determine active boundaries ---------------------------
tol = 1e-8;                         % Tolerance for when to consider a bound active
activeLowerX = lambdax > tol;       % Active lower x-bounds
activeUpperX = mux > tol;           % Active upper x-bounds
activeLowerY = lambday > tol;       % Active lower y-bounds
activeUpperY = muy > tol;           % Active upper y-bounds
activeLowerZ = lambdaz > tol;       % Active lower z-bounds
activeUpperZ = muz > tol;           % Active upper z-bounds

% ------------------------- Boundary points--------------------------------
markerSize = 10;
alpha = 1;

% Lower x
scatter3(solution(activeLowerX, 1), ...
    solution(activeLowerX, 2), ...
    solution(activeLowerX, 3), ...
    markerSize, "m", ...
    MarkerEdgeAlpha=alpha, ...
    MarkerFaceAlpha=alpha);

% Upper x
scatter3(solution(activeUpperX, 1), ...
    solution(activeUpperX, 2), ...
    solution(activeUpperX, 3), ...
    markerSize, "m", ...
    MarkerEdgeAlpha=alpha, ...
    MarkerFaceAlpha=alpha);

% Lower y
scatter3(solution(activeLowerY, 1), ...
    solution(activeLowerY, 2), ...
    solution(activeLowerY, 3), ...
    markerSize, "m", ...
    MarkerEdgeAlpha=alpha);

% Upper y
scatter3(solution(activeUpperY, 1), ...
    solution(activeUpperY, 2), ...
    solution(activeUpperY, 3), ...
    markerSize, "m", ...
    MarkerEdgeAlpha=alpha, ...
    MarkerFaceAlpha=alpha);

% Lower z
scatter3(solution(activeLowerZ, 1), ...
    solution(activeLowerZ, 2), ...
    solution(activeLowerZ, 3), ...
    markerSize, "m");

% Upper z
scatter3(solution(activeUpperZ, 1), ...
    solution(activeUpperZ, 2), ...
    solution(activeUpperZ, 3), ...
    markerSize, "m", ...
    MarkerEdgeAlpha=alpha, ...
    MarkerFaceAlpha=alpha);
% ----------------------- Additional plot settings ------------------------
axis equal;
xlabel("x");                % x-label
ylabel("y");                % y-label
zlabel("z");                % z-label
view(-70, 15);              % Horizontal - Vertical
legend("Originele mesh", "Vervormde mesh", "Actieve grenzen");  % Legend
title("Oplossing bol");     % Title

% Save the figure
if saveFigure, setFigParameters("deformation-sphere"); end
%% %%%%%%%%%%%%%%%%%%% Compare residual to refinement %%%%%%%%%%%%%%%%%%%%%
figure("Name", "Residuals of sphere", "NumberTitle", "off");
h = 0.07:0.01:0.20;                                 % Characteristic lenghts
hold on;
yscale("log");
for k = 1:length(h)
    history = load("Create_3D_meshes/Results/sphere_" + ...
        sprintf("%.2f", h(k)) + ".mat").history;    % Load mesh
    semilogy(history.residuals, "DisplayName", "h = " + h(k));  % Plot residual
end
xlabel("Iteratie")                                  % x-label
ylabel("Residu")                                    % y-label
title("Residuen")                                   % Title
legend()                                            % Legend
xl = xlim;                                          % Current limit
xlim([xl(1), xl(2) + 500]);                         % Extend the domain

% Save the figure
if saveFigure, setFigParameters("deformation-sphere_res"); end
%% %%%%%%%%%%%%%%%%%%%% Compare number of active bounds %%%%%%%%%%%%%%%%%%%
figure("Name", "Amount of active bounds, circle", "NumberTitle", "off");
hold on;

active_bounds = zeros(length(h), 1);
for k = 1:length(h)
    history = load("Create_3D_meshes/Results/sphere_" + ...
        sprintf("%.2f", h(k)) + ".mat").history;            % Load mesh

    tol = 1e-8;                                             % Tolerance to set active bounds
    active_lower = sum(history.lambda(:, end) > tol);       % Active lower bounds
    active_upper = sum(history.mu(:, end) > tol);           % Active upper bounds
    active_bounds(k) = active_lower + active_upper;         % Active bounds
end
scatter(h, active_bounds);                                  % Plot the results

xlabel("Karakteristieke lengte")                            % x-label
ylabel("Aantal actieve grenzen")                            % y-label
title("Vergelijking karakteristieke lengte, " + ...
    "aantal actieve grenzen")                               % Title

% Save figure
if saveFigure, setFigParameters("deformation-sphere_nrbounds"); end

%% %%%%%%%%%%%%%%%%%%% Compare with preconditioning %%%%%%%%%%%%%%%%%%%%%%%
% Load unpreconditioned results
spheres = {};
pspheres = {};
h = ["0.07", "0.10", "0.14", "0.20"]; % Characteristic lengths to compare

for k = 1:length(h)
    spheres{k} = load("Create_3D_meshes/Results/sphere_" + h(k) + ".mat");
    pspheres{k} = load("Create_3D_meshes/Results/sphere_" + h(k) + "_precond.mat");
end

figure("Name", "Residuals with preconditioning", "NumberTitle", "off");
hold on;
yscale("log");
c = get(gca, 'ColorOrder');                    % Get color palet
alpha = 0.7;                                   % amount of darkening

% Unpreconditioned
for k = 1:length(spheres)
    col = c(k, :);                              % Color
    semilogy(spheres{k}.history.residuals, "Color", col);     % Plot unpreconditioned residual
end

% Preconditioned
for k = 1:length(pspheres)
    col = c(k, :) * alpha;                      % Color
    semilogy(pspheres{k}.history.residuals, "Color", col);    % Plot preconditioned residual
end

legend("h = 0.07", "h = 0.10", "h = 0.14", "h = 0.20");       % Legend
xlabel("Iteratie");                             % x-label
ylabel("Residu");                               % y-label
title("Vergelijking residuen met en zonder preconditioning"); % Title

% Save results
if saveFigure, setFigParameters("3D-sphere");  end