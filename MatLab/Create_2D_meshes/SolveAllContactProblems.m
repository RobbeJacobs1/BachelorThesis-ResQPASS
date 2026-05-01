% This script solves all 2D contact problems
clc; clear all; close all;
%% %%%%%%%%%%%%%%%%%%%%%%%% Problems to solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("Create_2D_meshes/")

refinements = 0.05:0.01:0.2;                        % The characteristic lengts
test_files = strings(1, 3*length(refinements));       % Filenames
lb_x = zeros(size(test_files));
ub_x = zeros(size(test_files));
lb_y = zeros(size(test_files));
ub_y = zeros(size(test_files));
for k = 1:length(refinements)
    test_files(3*k-2) = "circle_" + sprintf('%.2f', refinements(k)); % Pick filename
    lb_x(3 * k - 2) = -0.8; % Lower x-bound
    ub_x(3 * k - 2) = 0.8;  % Upper x-bound
    lb_y(3 * k - 2) = -0.8; % Lower y-bound
    ub_y(3 * k - 2) = 0.8;  % Upper y-bound

    test_files(3*k-1) = "H_" + sprintf('%.2f', refinements(k));
    lb_x(3 * k - 1) = -0.4;
    ub_x(3 * k - 1) = 0.4;
    lb_y(3 * k - 1) = -1.5;
    ub_y(3 * k - 1) = 1.5;

    test_files(3*k) = "rectangle_" + sprintf('%.2f', refinements(k));
    lb_x(3 * k) = -1;
    ub_x(3 * k) = 1;
    lb_y(3 * k) = -1;
    ub_y(3 * k) = 1;
end
%%
% Solve every contact problem
for k = 1:length(test_files)
    test_files(k)
    solve_problem(test_files(k), lb_x(k), ub_x(k), lb_y(k), ub_y(k));
    solve_problem(test_files(k), lb_x(k), ub_x(k), lb_y(k), ub_y(k), true);
end