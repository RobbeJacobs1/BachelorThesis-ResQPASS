% This script aims to provide an example of the gradient in a minimum
% Determines whether the result should be saved.

% This corresponds to Example 2.2.2 (Figure 2.1, p. 5) in the accompanying paper.

% Decide wether the result should be formatted and saved.
saveFigure = false;

% Define function
f = @(x, y) (x.^2 + y.^2);

% Create mesh grid
n = 0.1;
[xGrid, yGrid] = meshgrid(-3:n:3, -2:n:2);

% Show the gradient
figure; hold on;
contourf(xGrid, yGrid, f(xGrid, yGrid), 30);
axis equal;

% Draw the boundary
theta = linspace(0, 2*pi, 200);
xc = cos(theta);
yc = sin(theta);
plot(xc, yc, 'w');
text(cos(pi/3) + 0.1, sin(pi/3) + 0.1, "$c(x)$", "Interpreter", "latex", ...
    "Color", "w")

% Plot optimal point
scatter(0, 0, 'co', "filled");
text(0.1, 0.1, "$x^*$", "Interpreter", "latex", "Color", "c")
xlabel("$x$", "Interpreter", "latex");
ylabel("$y$", "Interpreter", "latex");
title("$f(x, y) = x^2 + y^2$", "Interpreter", "latex");

if saveFigure setFigParameters("KKT voorbeeld interieur"); end