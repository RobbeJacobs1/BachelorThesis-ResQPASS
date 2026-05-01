function [x, history] = ResQPASS_altered(A,b,lb,ub, varargin)
% This implementation was obtained from https://github.com/AppliedMathUAntwerpen/ResQPASS/blob/main/Matlab/ResQPASSv3.m
% The only alterations made are to the way extra logging variables are
% stored and allowing the use of an initial value diffrent from 0.

% ResQPASS Residual Quadratic Programming Active Set Subspace algorithm.

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

% trueSol is used to check the error of the ResQPASS algorithm when the
% true solution is known. For more information, see the section on the
% history output. This is left empty by default.

% M1 is used as a preconditioner in the form of a function. Each
% iterations, M1(residual) will be calculated.

% maxInnerIt limits the amount of iterations used in the solver used for
% the subspace problem. This is by default set to 10.

% maxOuterIt limits the amount of outer iterations, meaning the amount of
% times the subspace problem is solved. This is by default min(m, n).

% doWarmStart is a logical variables which decides wether the algorithm
% should use warmstarting. This variable is true by default.

% doRecursiveQPAS is a logical variable as well which turns recursive 
% calculations in the QPAS subroutine on or off. This is true by default.

% resTol decides the tolerance that should be used to determine whether the
% algorithm has converged (when the norm of the residual goes below this
% variable). By default, resTol is 1e-8.

% cholTol gives the tolerance that should be used during the cholesky
% factorisation. This is by default set to 1e-10.

% displayResiduals can be used to print an output each iteration to monitor
% the pregression of the algorithm during longer calculations. This is
% false by default.

%%%%%%%%%%%%%%%%%%%%
% Output variables %
%%%%%%%%%%%%%%%%%%%%
% The first, and only necassary output variable is x. This variables 
% returns the solution ResQPASS delivered after it either converged or
% reached it's maximal number of iterations (maxOuterIt).

% history is a construct which contains all extra information about the
% iterations the algorithm made. This information will only be calculated
% if two input arguments are given.

% history.solutions contains a matrix with at column k, the solution at
% iteration k.

% history.objective contains the value of the objective function 
% ||A*x_k - b || at every iteration k. (here k is the k-th solution)

% history.residuals contains the norm of the residual at each iteration.

% history.errors is a vector containing the norm of the diffrence between
% the solution at iteration k and the actual solution provided in trueSol.

% history.lagrange_multipliers contains two variables (.lambda and .mu)
% which hold the lagrange multipliers at each iterations as columns of a
% matrix.

% history.innerIterations stores the amount of iterations the internal
% solver needed every iteration.

% history.recursiveError stores the recursion error for every inner QPAS 
% iteration.

% history.V contains the basismatrix of the Krylov-subspace at the time of
% convergence

% history.y contains the coordinates of the last solution in the Krylov 
% basis.

% history.outerIt is the amount of iterations the algorithm took until
% convergence

% history.mu and history.lambda are the Lagrange multipliers at the last
% iteration.

% history.converged is a logical variable which is true when the algorithm
% stopped because of convergence instead of reaching the maximal number of
% iterations.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default input handling %
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4               % Check for enough inputs
    error('Not enough input arguments.')
end

% Parse all arguments
[M,N] = size(A);            % Size of the problem
p = inputParser;
addParameter(p, "maxInnerIt", 10, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && mod(x,1)==0);
addParameter(p, "M1", [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, "maxOuterIt", [], ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && mod(x,1)==0);
addParameter(p, "doWarmStart", true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, "doRecursiveQPAS", true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, "resTol", 1e-8, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, "x0", [], ...
    @(x) isempty(x) || (isnumeric(x) && length(x) == size(A, 2) && isvector(x)));
addParameter(p, "cholTol", 1e-10, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, "trueSol", [], ...
    @(x) isempty(x) || (isnumeric(x) && length(x) == size(A, 2) && isvector(x)));
addParameter(p, "displayResiduals", false, @(x) islogical(x));

parse(p, varargin{:});
maxInnerIt = p.Results.maxInnerIt;
maxOuterIt = p.Results.maxOuterIt;
doWarmStart = p.Results.doWarmStart;
doRecursiveQPAS = p.Results.doRecursiveQPAS;
EPSres = p.Results.resTol;
EPSchol = p.Results.cholTol;
M1 = p.Results.M1;
x0 = p.Results.x0;
true_sol = p.Results.trueSol;
dispRes = p.Results.displayResiduals;

if isempty(maxOuterIt)      % Standard maxOuterIt
    maxOuterIt = min(M, N);
end

if isempty(lb)              % Create infinit lower bound if none provided
    lb = -Inf * ones(N, 1);
end

if isempty(ub)              % Create infinit upper bound if none provided
    ub = Inf * ones(N, 1);
end

doPreconditioning = ~isempty(M1);   % Check for preconditioning 

if ~isempty(x0)
    if size(x0(lb > 0), 1) > 0 || size(x0(ub < x0), 1) > 0      % Check if x0 is feasible
        warning("The given x0 is not feasible, the algorithm will continue with " + ...
            "a diffrent feasible x0");
        x0 = [];
    end
end

if isempty(x0)                              % Create x0 in the middle of both bounds
    lb_inf = (lb == -Inf & ub < Inf)';
    ub_inf = (lb > -Inf & ub == Inf)';
    no_inf = (lb > -Inf & ub < Inf)';
    
    x0 = zeros(N, 1);
    x0(lb_inf) = ub(lb_inf) - 1;
    x0(ub_inf) = lb(ub_inf) + 1;
    x0(no_inf) = (ub(no_inf) + lb(no_inf))/2;
end

if nargout > 1
    logging = true;
    history.solutions = zeros(N, maxOuterIt);
    history.residuals = zeros(maxOuterIt, 1);
    history.objective = zeros(maxOuterIt, 1);
    if ~isempty(true_sol)
        history.errors = zeros(maxIter, 1);
    end
    history.lagrange_multipliers.lambda = zeros(N, maxOuterIt);
    history.lagrange_multipliers.mu = zeros(N, maxOuterIt);
    history.innerIterations = zeros(maxOuterIt, 1);
    history.recursiveError = [];
end

converged = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translate the problem for x0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translate the problem using x = z + x0 and solving for z starting from 
% z = 0, which is feasible.
b = b - A * x0;
lb = lb - x0;
ub = ub - x0;

%%%%%%%%%%%%%%%%%%
% Initialisation %
%%%%%%%%%%%%%%%%%%
v = A'*b;                   % Initial residual
res = norm(v);           % Norm of residual
V = v/res;               % Initial basis
AV = A*V;                   % Helper variable
L = sqrt(AV'*AV);           % Lower cholesky factor of hessian V'A'AV (=LL')
f = -b'*AV;                 % Linear factor
workingSet = []; y = [];    % Working set for QPAS
lu = [-lb,ub];              % Right hand side for qpas inequality
VV = [-V;V];                % Inequality matrix for qpas TODO: REMOVE THE USE OF V SOMEHOW?

%%%%%%%%%%%%%%
% Outer loop %
%%%%%%%%%%%%%%
for outerIt = 1:maxOuterIt
    if dispRes, disp("Iteration: " + outerIt + ", Residual: " + res); end   % Display residuals
    % Bidiag
    U(:,outerIt) = AV(:,end);
    for j = 1:outerIt-1
        B(j,outerIt) = dot(U(:,j),U(:,end));
        U(:,end) = U(:,end) - B(j,outerIt)*U(:,j);
    end
    B(outerIt, outerIt) = norm(U(:,end));
    U(:,end) = U(:,end)/B(outerIt, outerIt);
    
    % diff{outerIt} = L-B';

    L = B';

    %%%%%%%%%%%%%%%%%%%%
    % Inner loop (QPAS)%
    %%%%%%%%%%%%%%%%%%%%
    if ~doWarmStart         % Reset working set  and y when cold-starting
        workingSet = [];
        y = zeros(outerIt - 1, 1);
    end
    if ~logging             % Standard QPAS
        [y,workingSet,lagMultActive] = ...
            qpasCholeskyv2(L,f,VV,lu,[y;0],workingSet,maxInnerIt,doRecursiveQPAS);
    else                    % QPAS with computation of recursion error
        [y,workingSet,lagMultActive,innerIters,recursiveErrorIter] = ...
            qpasCholeskyv2(L,f,VV,lu,[y;0],workingSet,maxInnerIt,doRecursiveQPAS);
        history.recursiveError = [history.recursiveError, recursiveErrorIter];
    end
    % Lagrange multipliers
    lagMult = zeros(2*N,1);
    lagMult(workingSet) = lagMultActive;
    lambda = lagMult(1:N);          % Lower bounds
    mu = lagMult(N+1:end);          % Upper bounds

    %%%%%%%%%%%%
    % Residual %
    %%%%%%%%%%%%
    temp = (AV*y-b);                % Used for objective and residual
    v = A'*temp - lambda + mu;      % Calculate new residual
    if doPreconditioning            % Apply preconditioner
        v = M1(v);
    end
    res = norm(v);                  % Calculate norm of residual

    %%%%%%%%%%%%%%%%%%%%%
    % Log all variables %
    %%%%%%%%%%%%%%%%%%%%%
    if logging
        x = V * y + x0;
        history.solutions(:, outerIt) = x;
        history.residuals(outerIt) = res;
        history.objective(outerIt) = temp' * temp;
        if ~isempty(true_sol)
            history.errors(outerIt) = norm(true_sol - x);
        end
        history.lagrange_multipliers.lambda(:, outerIt) = lambda;
        history.lagrange_multipliers.mu(:, outerIt) = mu;
        history.innerIterations(outerIt) = innerIters;
    end


    %%%%%%%%%%%%
    % Stopping %
    %%%%%%%%%%%%
    if res < EPSres        % Small residual
        converged = true;
        break
    elseif L(end,end) < EPSchol     % Breakdown of choleksy (happy breakdown)
        converged = true;
        break
    end
    
    %%%%%%%%%%%
    % Updates %
    %%%%%%%%%%%
    V = [V v/res];     % Basis expansion
    VV = [VV, [-V(:,end);V(:,end)]];                % Inequality matrix for qpas TODO: REMOVE THE USE OF V SOMEHOW?
    AV = [AV,A*V(:,end)]; %TODO: we can do it without (through L only)?
    % U12 = L\(AV(:,1:outerIt)'*AV(:,outerIt+1));
    
    % Cholesky
    % L = [L, zeros(outerIt,1);
    %         U12',sqrt(AV(:,outerIt+1)'*AV(:,outerIt+1)-U12'*U12)];
    

    % % Bidiag
    % U(:,outerIt) = AV(:,end);
    % for j = 1:outerIt-1
    %     B(j,outerIt) = dot(U(:,j),U(:,end));
    %     U(:,end) = U(:,end) - B(j,k)*U(:,j);
    % end
    % B(outerIt, outerIt) = norm(U(:,end));
    % U(:,end) = U(:,end)/B(outerIt, outerIt);
    % L = B;

    f = [f;-b'*AV(:,end)];      % Linear term
end

x = V * y + x0;         % Find the untranslated solution
if logging
    history.solutions = history.solutions(:, 1:outerIt);
    history.residuals = history.residuals(1:outerIt);
    history.objective = history.objective(1:outerIt);
    if ~isempty(true_sol)
        history.errors = history.errors(1:outerIt);
    end
    history.lagrange_multipliers.lambda = history.lagrange_multipliers.lambda(:, 1:outerIt);
    history.lagrange_multipliers.mu = history.lagrange_multipliers.mu(:, 1:outerIt);
    history.innerIterations = history.innerIterations(1:outerIt);
    history.V = V;
    history.y = y;
    history.outerIt = outerIt;
    history.mu = mu;
    history.lambda = lambda;
    history.converged = converged;
end
if ~converged
    warning("ResQPASS did not converge in %i iterations", maxOuterIt)
end
end