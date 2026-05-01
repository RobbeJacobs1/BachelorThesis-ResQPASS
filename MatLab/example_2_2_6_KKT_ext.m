% This script aims to give an example when the optimal point lies on the
% edge of the boundary. This is then compared to the gradient of the
% boundary and the gradient of the objective function.

% This corresponds to Example 2.2.6 (Figure 2.2, p. 10) in the accompanying paper.

% Determine whether the figure should be saved.
saveFigure = false;

% Define the function
f = @(x, y) 0.5 * (x.^2 + y.^2);

% Define the meshgrid
n = 0.1;
[xGrid, yGrid] = meshgrid(-2:n:4, -2:n:2);

% Make a contourplot of the objective function
figure; hold on;
contourf(xGrid, yGrid, f(xGrid, yGrid), 30);
axis equal;

% Plot the boundary
theta = linspace(0, 2*pi, 200);
xc = 2 + cos(theta);
yc = sin(theta);
plot(xc, yc, 'w', 'LineWidth', 2);

% Calculate the gradient
x0 = 1; y0 = 0;
grad = [x0, y0];

% Plot the gradient of the objective function
quiver(x0, y0, grad(1), -grad(2), norm(grad), 'r');

% Determine the gradient of the boundary function
tangent_dir = [grad(2), -grad(1)];

% Plot the gradient of the boundary function
t = linspace(-1, 1, 100);
xt = x0 + t * tangent_dir(1);
yt = y0 + t * tangent_dir(2);
xline(xt, "m", "LineWidth", 1);


% Plot optimal point
scatter(x0, y0, 'co', "filled");
text(2 + cos(pi/3), sin(pi/3) + 0.1, '$c(x)$', "Interpreter", "latex", ...
    'Color', 'w');
text(x0+0.1, y0+0.1, "$x^*$", "Interpreter", "latex", "Color", "c")

% Gradient labels
text(x0+1, y0+0.1, "$\nabla f(x)$", "Interpreter", "latex", ...
    "Color", "r");
text(1.1, 1.5, "$T_{x^*}$", "Interpreter", "latex", ...
    "Color", "m")

% Minimal point
scatter(0, 0, "go", "filled")
text(0.1, 0.1, "$x_{min}$", "Interpreter", "latex", ...
    "Color", "g")
title("$f(x, y) = x^2 + y^2$", "Interpreter", "latex")
xlabel("$x$", "Interpreter", "latex");
ylabel("$y$", "Interpreter","latex");

% Save the figure
if saveFigure setFigParameters("KKT voorbeeld"); end