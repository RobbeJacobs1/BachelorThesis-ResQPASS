% This script solve all given contact problems. 
close all; clear all; clc;
%% %%%%%%%%%%%%%%%%%%%%%%%% Problems to solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("Create_3D_meshes/")

refinements = flip(0.07:0.01:0.20);             % Flip to start with eassiest problems for fast results
test_files = strings(1, length(refinements));

for k = 1:length(refinements)
    test_files(k) = "sphere_" + sprintf('%.2f', refinements(k));
end

%%
% Determine the bounds for the problems
lb_x = -0.8 * ones(size(test_files));
ub_x = 0.8 * ones(size(test_files));
lb_y = -0.8 * ones(size(test_files));
ub_y = 0.8 * ones(size(test_files));
lb_z = -0.8 * ones(size(test_files));
ub_z = 0.8 * ones(size(test_files));

for k = 1:length(test_files)
    test_files(k)
    % Solve the problem without preconditioner
    solve_problem(test_files(k), ...
        lb_x(k), ub_x(k), ...
        lb_y(k), ub_y(k), ...
        lb_z(k), ub_z(k), false);

    % Solve the problem with preconditioner
    solve_problem(test_files(k), ...
        lb_x(k), ub_x(k), ...
        lb_y(k), ub_y(k), ...
        lb_z(k), ub_z(k), true);
end
