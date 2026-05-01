function [xk, workingSet, lagMult, innerIt, errorRecursive] = qpasCholeskyv2(L,f,C,d,x0,workingSet,maxiter, ...
    recursiveFlag,showWarnings)
% This code was obtained from https://github.com/AppliedMathUAntwerpen/ResQPASS/blob/main/Matlab/qpasCholeskyv2.m

%qpasCholesky quadratic programming active-set with lower-Cholesky input.
%Subroutine of the ResQPASS method.
%
%   xk = qpasCholesky(L,f,C,d,x0) attempts to solve the quadratic
%   programming problem: 
%       min 0.5*x'*L*L'*x + f'*x subject to: C*x <= d
%        x
%   where L is the lower-Cholesky factor of the Hessian. The problem is 
%   solved with an active-set algorithm (qpas) with initial starting point x0. 
%   x0 should be feasible.
%
%   xk = qpasCholesky(L,f,D,d,x0,workingSet) warm-starts the problem with
%   some initial working set. (default: [])
%
%   xk = qpasCholesky(L,f,C,d,x0,workingSet,maxiter) limits the number of
%   iterations to maxiter. (default: 10*(numberOfVariables + numberOfConstraints))
%
%   xk = qpasCholesky(L,f,C,d,x0,workingSet,maxiter,recursiveFlag) uses a
%   recursion relation to compute C*xk if the flag is set to true. This is
%   unstable but faster. (default: false)
%
%   xk = qpasCholesky(L,f,C,d,x0,workingSet,maxiter,recursiveFlag,showWarnings)
%   shows warnings if set to true. (default: false)
%   
%   [xk, workingSet] = qpasCholesky(___) also returns the final working
%   set for use in warm-starting.
%
%   [xk, workingSet, lagMult] = qpasCholesky(___) returns the lagrange
%   multipliers of the active constraints, corresponding to the indices in
%   the working set.
%
%   [xk, workingSet, lagMult, innerIt] = qpasCholesky(___) returns the number
%   of iterations performed
%
%   [xk, workingSet, lagMult, innerIt, errorRecursive] = qpasCholesky(___) 
%   returns the recursion error of C*xk for every iteration.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default input handling %
%%%%%%%%%%%%%%%%%%%%%%%%%%
EPSzero = 1e-9;        % Default tolerance

if nargin < 5
    error('Not enough input arguments.')
elseif nargin == 5
    workingSet = [];
    maxiter = 5; % Default of MATLAB's quadprog (active-set)
    recursiveFlag = false;
    showWarnings = false;
elseif nargin == 6
    maxiter = 5; % Default of MATLAB's quadprog (active-set)
    recursiveFlag = false;
    showWarnings = false;
elseif nargin == 7
    recursiveFlag = false;
    showWarnings = false;
elseif nargin == 8
    showWarnings = false;
elseif nargin > 9
    error('Too many input arguments.')
end

%%%%%%%%%%%%%%%%%%
% Initialisation %
%%%%%%%%%%%%%%%%%%
m = size(C,1);                  % Number of constraints
Ginv_f = (L')\(L\f);            % Helper G^{-1}c 
workingSetCompl = MY_setdiff(1:m,workingSet);  %Complement of the working set
xk = x0;                        % Initial guess
C_xk = C*xk;                    % Helper
errorRecursive = [];            % Error in recursion of C*xk

% Initialisation for the first iteration
pk = -(xk + Ginv_f);        % Search direction (without Lagrange multipliers)
if isempty(workingSet)      %Cold start
    lagMult = [];           % Lagrange multipliers (none for empty working set)
    Linv_Ckt = [];          % Helper L^{-1}C_{W_k}
    Ginv_Ckt = [];          % Helper G^{-1}C_{W_k}
    Q = []; R = [];         % QR factors
    rhs = [];               % Right-hand side of the QPAS system
else                        % Warm start
    Linv_Ckt = linsolve(L,(C(workingSet,:)'), struct('LT', true));
    Ginv_Ckt = linsolve(L', Linv_Ckt, struct('UT', true));            % Helper G^{-1}C_{W_k}^T
%     schur = C(workingSet,:)*Ginv_Ckt;   % Schur Complement
    [Q,R] = qr(Linv_Ckt);                  % QR factors
    rhs = C(workingSet,:)*(-pk);        % Right-hand side of the QPAS system
        
    lagMult = -linsolve(R, ...
        linsolve(R',rhs, struct('LT',true)), ...
        struct('UT', true));     % Lagrange multipliers
    pk = pk - Ginv_Ckt*lagMult;         % Search direction (with Lagrange multipliers)
end

innerIt = 1;
invalidStopFlag = true;
while true
    %%%%%%%%%%%%%%%%%
    % Step possible %
    %%%%%%%%%%%%%%%%%
    if norm(pk) > EPSzero
        invalidStopFlag = true;
        %   Find & add the blocking contraint to the working set
        %   Find indices of all constraints not in the working set, where a_i^T * p > 0 
        C_pk = C*pk;
        C_pk_pos = find(C_pk > 0);
        indices  = MY_intersect(workingSetCompl,C_pk_pos);
        
        alphas = C_xk(indices);
        alphas = d(indices)-alphas;
        alphas = alphas./C_pk(indices);           % All possible stepsizes
        [alpha,min_alpha_idx] = min(alphas);    % Smallest stepsize
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % Blocking constraint %
        %%%%%%%%%%%%%%%%%%%%%%%
        if alpha<=1                             % Blocking constraint
            xk = xk + alpha*pk;                 % Update

            bcidx = indices(min_alpha_idx);     % Blocking index
            workingSet = [workingSet,bcidx];    % Add blocking index to working set
            workingSetCompl(workingSetCompl == bcidx) = [];     %Remove the blocking index from the complement

            % Update helpers, QR-factorization and rhs of qpas-system
            Linv_Ckt = [Linv_Ckt, linsolve(L,(C(bcidx,:)'), struct('LT', true))];
            Ginv_Ckt = [Ginv_Ckt, linsolve(L',Linv_Ckt(:,end), struct('UT', true))];
%             column = C(workingSet(1:end-1), :)*Ginv_Ckt(:,end);
%             row = C(bcidx, :)*Ginv_Ckt;
            [Q,R] = qrinsert(Q,R,size(R,2)+1,Linv_Ckt(:,end),'col');
%             [Q,R] = qrinsert(Q,R,size(Q,1)+1,row,'row');
            rhs = [rhs; C(bcidx,:)*(xk + Ginv_f)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % No blocking constraint %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            alpha = 1;          % Used in recursive computation of C*xk
            xk = xk + pk;
        end

        % Compute the new helper C*xk
        if recursiveFlag    % Recursively
            C_xk = C_xk + alpha*C_pk;
            if nargout > 4
                errorRecursive(innerIt)=norm(C_xk-C*xk);     % Error because of recursion
            end
        else                % Directly
            C_xk = C*xk;    
        end
    
    %%%%%%%%%%%
    % No step %
    %%%%%%%%%%%
    else
        invalidStopFlag = false;
        [minLagMult,minLagMult_idx] = min(lagMult); % Smallest lagrange multiplier

        %%%%%%%%%%%%%%%%%%%%%
        % Remove constraint %
        %%%%%%%%%%%%%%%%%%%%%
        if minLagMult < 0
            if innerIt >= maxiter
                if showWarnings
                    % Show a warning when the method did not converge within maxiter iterations
                    warning('QPAS did not converge in %i iterations', maxiter)
                end
                return
            end
            % Move the constraint to the non-working set
            workingSetCompl = [workingSetCompl, workingSet(minLagMult_idx)]; 
            workingSet(minLagMult_idx) = [];
            
            % Delete from helper, QR-factorization and rhs of qpas-system  
            
            Linv_Ckt(:,minLagMult_idx) = [];     
            Ginv_Ckt(:,minLagMult_idx) = [];     
            [Q,R] = qrdelete(Q,R,minLagMult_idx,'col');
%             [Q,R] = qrdelete(Q,R,minLagMult_idx,'row');
            rhs(minLagMult_idx) = [];

        %%%%%%%%%%%%%%%%%%
        % Solution found %
        %%%%%%%%%%%%%%%%%%
        else
            return
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % New search direction %
    %%%%%%%%%%%%%%%%%%%%%%%%
    pk = -(xk + Ginv_f);        % Search direction (without Lagrange multipliers)
    if isempty(workingSet)
        lagMult = [];           % Lagrange multipliers (none for empty working set)
        Linv_Ckt = [];          % Helper L^{-1}C_{W_k}
        Ginv_Ckt = [];          % Helper G^{-1}C_{W_k}
        Q = []; R = [];         % QR factors
        rhs = [];               % Right-hand side of the QPAS system
    else      
%         lagMult = -linsolve(R,(Q'*rhs), ...
%                     struct('UT', true));        % Lagrange multipliers
            lagMult = -linsolve(R, ...
        linsolve(R',rhs, struct('LT',true)), ...
        struct('UT', true));     % Lagrange multipliers
        pk = pk - Ginv_Ckt*lagMult;             % Search direction (with Lagrange multipliers)
    end

innerIt = innerIt+1;

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast intersection and difference for sets %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%From https://nl.mathworks.com/matlabcentral/answers/53796-speed-up-intersect-setdiff-functions

function C = MY_intersect(A,B)
 P = zeros(1, max(max(A),max(B)) ) ;
 P(A) = 1;
 C = B(logical(P(B)));
end

function Z = MY_setdiff(X,Y)
  check = false(1, max(max(X), max(Y)));
  check(X) = true;
  check(Y) = false;
  Z = X(check(X));  
end