function [sol, history] = ResQPASS_quadprog(A, b, lb, ub, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required input variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The only required input variables are does which specify the problem
% solved, which is min  || Ax - b ||.
%                lb<x<ub
% Where A is a real valued m x n matrix, b a vector in R^m and the bounds
% of the problem lb, ub in R^n.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional input variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 determines the first guess the algorithm starts from. By default, this
% is the zero-vector.

% sol is used to check the error of the ResQPASS algorithm when the
% true solution is known. For more information, see the section on the
% history output. This is left empty by default.

% preConditioner is used as a preconditioner in the form of a function.
% Each iterations, this function will be applied to the residual.

% maxIter limits the amount of outer iterations, meaning the amount of
% times the subspace problem is solved. This is by default max(m, n).

% warmStart is a logical variables which decides whether the algorithm
% should use warmstarting. This variable is true by default.

% tol decides the tolerance that should be used to determine whether the
% algorithm has converged (when the norm of the residual goes below this
% variable). By default, tol is 1e-8.

% algorithm is a string corresponding to the different algorithms quadprog
% can utilise. This is either active-set or interior-point (active-set by
% default). Note that interior-point cannot utilise initial guesses.

%%%%%%%%%%%%%%%%%%%%
% Output variables %
%%%%%%%%%%%%%%%%%%%%
% The first, and only necassary output variable is sol. This variables 
% returns the solution ResQPASS delivered after it either converged or
% reached it's maximal number of iterations (maxIter).

% history is a construct which contains all extra information about the
% iterations the algorithm made. This information will only be calculated
% if two input arguments are given.

% history.solutions contains a matrix with at column k, the solution at
% iteration k.

% history.objective contains the value of the objective function 
% ||A*x_k - b || at every iteration k. (here k is the k-th solution)

% history.residuals contains the norm of the residual at each iteration.

% history.errors is a vector containing the norm of the diffrence between
% the solution at iteration k and the actual solution provided in sol.

% history.lambda contains a matrix with at each column k, the lagrange
% multiplier lambda for the lower bound at iteration k.

% history.mu contains a matrix with at each column k, the lagrange
% multiplier mu for the upper bound at iteration k.

% history.innerIterations stores the amount of iterations the internal
% solver needed every iteration

% history.V contains the basismatrix of the Krylov-subspace at the time of
% convergence

% history.y contains the coordinates of the last solution in the Krylov 
% basis.

% history.iters is the amount of iterations the algorithm took until
% convergence

% history.converged is a logical variable which is true when the algorithm
% stopped because of convergence instead of reaching the maximal number of
% iterations.

    %%%%%%%%%%%%%%%%%%%%%% Parse input parameters %%%%%%%%%%%%%%%%%%%%%%%%%
    [m, n] = size(A);

    doHistory = nargout > 1;      % If history is passed, extra intermediate information is saved.

    p = inputParser;
    addParameter(p, "maxIter", [], @(x) isnumeric(x) && isscalar(x) && x > 0);   % Maximal amount of outer iterations
    addParameter(p, "tol", 1e-8, @(x) isnumeric(x) && isscalar(x) && x > 0);     % When KKT-condition go below tol, stop
    addParameter(p, "sol", [], ...
        @(x) isempty(x) || (isequal(size(x), [n, 1]) && isnumeric(x)));          % Solution to compare the used solution to
    addParameter(p, "preConditioner", [], ...
        @(x) isempty(x) || isa(x, "function_handle"));                           % Preconditioner as a funtion handle
    addParameter(p, "x0", [], @(x) isempty(x) || ...
        (isnumeric(x) && isequal(size(x), [n, 1])));                             % Add initial guess
    addParameter(p, "warmStart", false, @(x) islogical(x));
    addParameter(p, "algorithm", "active-set");

    parse(p, varargin{:});
    maxIter = p.Results.maxIter;
    tol = p.Results.tol;
    true_sol = p.Results.sol;
    M =  p.Results.preConditioner;
    x0 = p.Results.x0;
    warmStart = p.Results.warmStart;
    algorithm = p.Results.algorithm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%
    % If no maximum amount of iterations is provided, set it to the largest
    % dimension of A.
    if isempty(maxIter)
        maxIter = max(m, n);
    end

    % If one of the bounds is empty, set it to infinity.
    if isempty(lb)
        lb = -Inf * ones(n, 1);
    end

    if isempty(ub)
        ub = Inf * ones(n, 1);
    end

    % Check if x0 is feasible, if it is not, an other x0 will be used
    if ~isempty(x0)
        if size(x0(lb > x0), 1) > 0 || size(x0(ub < x0), 1) > 0
            warning("The given x0 is not feasible.")
            x0 = [];
        end
    end

    % Pick x0 between the bounds, if one of the bounds is infinite, x0 is
    % offset away from the bound. If both bounds are infinite, x0 is set to
    % 0.
    if isempty(x0)
        lb_inf = (lb == -Inf & ub < Inf)';
        ub_inf = (lb > -Inf & ub == Inf)';
        no_inf = (lb > -Inf & ub < Inf)';
        
        x0 = zeros(n, 1);
        x0(lb_inf) = ub(lb_inf) - 1;
        x0(ub_inf) = lb(ub_inf) + 1;
        x0(no_inf) = (ub(no_inf) + lb(no_inf))/2;
    end
    
    % Apply a linear transformation to make the initial guess be 0
    b = b - A * x0;
    lb = lb - x0;
    ub = ub - x0;

    % Find where the bounds matter
    finite_lb = isfinite(lb);
    finite_ub = isfinite(ub);
    nr_active_lb = sum(finite_lb);
    nr_active_ub = sum(finite_ub);

    r = A' * b;                     % Compute initial residual

    if ~isempty(M)                  % Apply pre-conditioner
        r = M(r);
    end
    norm_r = norm(r);
    V = zeros(n, maxIter + 1);      % Initialise the Krylov basis
    V(:, 1) = r/norm_r; 
    y = [];                         % The solution in the basis will be stored in y

    fval = -A' * b;                 % Helper variable for quadprog
    bbound = [-lb(finite_lb); ub(finite_ub)];             % Bound in the form used for quadprog

    if doHistory                                % Initialise all history fields for logging
        history.solutions = zeros(n, maxIter);
        history.residuals = zeros(maxIter, 1);
        history.objective = zeros(maxIter, 1);
        if ~isempty(true_sol)
            history.errors = zeros(maxIter, 1);
        end
        history.lambda = zeros(n, maxIter);
        history.mu = zeros(n, maxIter);
        history.converged = false;
        history.innerIterations = zeros(maxIter, 1);
    end

    quadprog_options = optimoptions('quadprog', 'Display', 'off', ...   % Turn off quadprog messages
                                    "Algorithm", algorithm);         % Use active set
    

    % Initialse lagrange multipliers
    lambda = zeros(n, 1);
    mu = zeros(n, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter = 1:maxIter
        % ------------ Solve the problem in the Krylov subspace -----------
        % Put the problem into a form quadprog can use
        AV = A * V(:, 1:iter);
        H = AV' * AV;
        f = V(:, 1:iter)' * fval;
        Abound = [-V(finite_lb, 1:iter); V(finite_ub, 1:iter)];
        if warmStart
            y = [y; 0];
        else
            y = zeros(size(y, 1) + 1, 1);
        end

        [y, ~, ~, output, lagrange_multipliers] = quadprog(H, f, Abound, bbound, ...
            [], [], [], [], y, quadprog_options);
   
        % ---------------- Determine lagrange multipliers -----------------
        current_idx = 1;
        if nr_active_lb ~= 0
            lambda(finite_lb) = lagrange_multipliers.ineqlin(1:nr_active_lb);   % Lower bound
            current_idx = nr_active_lb + 1;
        end
        if nr_active_ub ~= 0
            mu(finite_ub) = lagrange_multipliers.ineqlin(current_idx:end);
        end

        % ------------------------ Expand the basis -----------------------
        sol = V(:, 1:iter) * y;                 % Compute the solution (still translated)
        r = A' * (A * sol - b) - lambda + mu;   % Compute next residual
        if ~isempty(M)                          % Apply preconditioner
            r = M(r);
        end
        norm_r = norm(r);

        V(:, iter + 1) = r / norm_r;            % Expand the basis

        % ------------------ Store everything in history ------------------
        if doHistory
            untranslated_sol = sol + x0;        % Translate the solution back
            history.solutions(:, iter) = untranslated_sol;   
            history.residuals(iter) = norm_r;
            history.objective(iter) = norm(A*sol - b)^2;
            if ~isempty(true_sol)
                history.errors(iter) = norm(untranslated_sol - true_sol);
            end
            history.lambda(:, iter) = lambda;
            history.mu(:, iter) = mu;
            history.innerIterations(iter) = output.iterations;
        end

        % ----------------------- Check convergence -----------------------
        if norm_r < tol
            history.converged = true;
            break
        end
    end

    %%%%%%%%%%%% Adjust history values for amount of iterations %%%%%%%%%%%
    if doHistory 
        if history.converged
            sol = untranslated_sol;
            history.iters = iter;
            history.solutions = history.solutions(:, 1:iter);
            history.residuals = history.residuals(1:iter);
            history.objective = history.objective(1:iter);
            if ~isempty(true_sol)
                history.errors = history.errors(1:iter);
            end
            history.mu = history.mu(:, 1:iter);
            history.lambda = history.lambda(:, 1:iter);
            history.V = V(:, 1:iter + 1);
            history.innerIterations = history.innerIterations(1:iter);
            history.y = y;
            return
        end
    end

    sol = sol + x0;
    if iter == maxIter
        warning("ResQPASS did not converge");
    end
    
end