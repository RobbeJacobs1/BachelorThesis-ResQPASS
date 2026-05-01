function [x, history] = lsqr_self(A, b, x0, maxit, tol, M)
    % This is a simple implementation of the LSQR-algoritm where the
    % problem is solved in the subspaces using the built-in lsqr-function
    % of MatLab. This is mainly used to get acces to the intermitent
    % iterations.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Required input variables %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The only required input variables are A and b. These are used to
    % define the problem that will be solved:
    %       min ||Ax - b||

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optional input variables %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x0 represents the initial value. It is 0 by default

    % maxit limits the amount of iterations are performed, the number of
    % columns of A by default.

    % tol gives the tollerance to determine when the algorithm should
    % terminate (when the norm of the residual goes below it). 
    % This is 1e-8 by default.

    % M represents a preconditioner. This is a function which will be
    % applied to the residual at each iteration. This is empty by default.

    %%%%%%%%%%%%%%%%%%%%
    % Output variables %
    %%%%%%%%%%%%%%%%%%%%
    % There are two output variables: x and history.
    % x stores the solution at the last iteration.

    % history is a construct containting information about the iterations
    % performed.

    % history.residuals is a vector containting the norm of the residuals
    % at each iteration.

    % history.solutions is a matrix with as each column k, the solution at
    % iteration k.

    % history.V contains the basis constructed at each iterations
    % (normalized residuals) as columns.

    % history.objective keeps the value of the objective function at each
    % iteration. This value is ||A x_k - b|| with x_k the solution at
    % iteration k.

    % history.iters gives the number of iterations needed before
    % convergence.

    % history.converged is a logic variable which is true if the algorithm
    % terminate because of convergence, false if the maximal number of
    % iterations was reached.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Handle input variables %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    n = size(A, 2);             % Size of A
    if nargin < 3               % Initial guess
        x0 = zeros(n, 1);
    end
    if nargin < 4 % Max iterations
        maxit = n;
    end
    if nargin <5 % Tolerance
        tol = 1e-8;
    end
    if nargin < 6 % Preconditioner
        M = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize variables %
    %%%%%%%%%%%%%%%%%%%%%%%%    
    r = A' * b; % First residual
    if ~isempty(M) % Apply preconditioner
        r = M(r);
    end

    V = zeros(n, maxit + 1); % Create basis matrix
    V(:, 1) = r / norm(r); % Append the first residual vector
    y = [];
    % Initialize history object
    history.residuals = zeros(maxit, 1);
    history.solutions = zeros(n, maxit);
    history.objective = zeros(maxit, 1);
    history.converged = false;
    for k = 1:maxit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve the problem in the subspace %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y = [y; 0]; % New initial guess for warm starting
        [y, ~] = lsqr(A * V(:, 1:k), b - A*x0, tol, [], [], [], y); % Solve the subspace problem
        x = x0 + V(:, 1:k)*y; % Calculate the solution corresponding to the basis
        v = A' * (A * x - b); % Calculate next residual
        if ~isempty(M)
            v = M(v); % Apply preconditioner
        end

        for j = 1:k
            v = v - (V(:,j)' * v) * V(:,j); % Reorthogonalisation
        end

        V(:, k + 1) = v/norm(v); % Append the residual vector

        %%%%%%%%%%%%%%%%%%%%%%%%
        % Keep track of values %
        %%%%%%%%%%%%%%%%%%%%%%%%
        history.residuals(k) = norm(v); % Residual
        history.solutions(:, k) = x; % Solutions
        history.objective(k) = norm(A*x - b)^2; % Objective
        
        %%%%%%%%%%%%%%%%%%%%%
        % Check convergence %
        %%%%%%%%%%%%%%%%%%%%%
        if norm(v) < tol
            history.converged = true;
            break;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Shorten every variable used for logging to number of iterations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    history.residuals = history.residuals(1:k);
    history.V = V(:, 1:k+1);
    history.solutions = history.solutions(:, 1:k);
    history.objective = history.objective(1:k);
    history.iters = k;
end