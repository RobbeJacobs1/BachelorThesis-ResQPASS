% This function solves a specific contact problem
function [] = solve_problem(filename, lb_x, ub_x, lb_y, ub_y, preconditioning)
    % Initialise values
    if nargin < 6
        preconditioning = false;
    end

    filepath = "Create_2D_meshes/Parsed_meshes/" + filename + ".mat";
    
    % Load files
    load(filepath, "nodes", "triangles");

    %% Construct graph Laplacian/stiffness matrix
    nV = size(nodes,1);                  % number of vertices
    
    % Collect all edges from triangles
    E = [triangles(:,[1 2]);             % edges between vertices 1-2
         triangles(:,[2 3]);             % edges between vertices 2-3
         triangles(:,[3 1])];            % edges between vertices 3-1
    
    % Sort edges so (i,j) and (j,i) are treated identically
    E = sort(E, 2);
    
    % Remove duplicate edges
    E = unique(E, 'rows');
    
    d = 2;                                 % spatial dimension
    Ndof = nV*d;                           % total DOFs
    
    Q = zeros(Ndof, Ndof);
    k0 = 1;                                % Spring constant 
    
    nr_edges = size(E, 1);
    
    % Loop over all edges to build stiffness matrix
    for e = 1:nr_edges
        i = E(e,1);                     % first vertex of edge
        j = E(e,2);                     % second vertex of edge
        
        % Calculate direction vector and length
        r = nodes(j,:) - nodes(i,:);          % vector from i to j
        l = norm(r);                    % edge length
        n_ij = r / l;                   % unit direction vector
        
        % Local stiffness matrix (rank-1 matrix for spring)
        Ke = k0 * (n_ij' * n_ij); % Factor k0 / l makes it move down (k0 has no impact)
        
        % DOF indices for vertices i and j
        Ii = (i-1)*d + (1:d);           % DOF indices for vertex i [2i-1, 2i]
        Ij = (j-1)*d + (1:d);           % DOF indices for vertex j [2j-1, 2j]
        
        % Assemble into global Laplacian matrix
        Q(Ii,Ii) = Q(Ii,Ii) + Ke;       % diagonal block for vertex i
        Q(Ij,Ij) = Q(Ij,Ij) + Ke;       % diagonal block for vertex j
        Q(Ii,Ij) = Q(Ii,Ij) - Ke;       % off-diagonal coupling i-j
        Q(Ij,Ii) = Q(Ij,Ii) - Ke;       % off-diagonal coupling j-i
    end

    %% Solve optimization problem with box constraints    
    % Convert the bounds to differences 
    lb = [lb_x * ones(nV, 1) lb_y * ones(nV, 1)] - nodes;
    ub = [ub_x* ones(nV, 1) ub_y* ones(nV, 1)] - nodes;
    
    % Reshape to get [x y; ...] --> [x ; y ; ...]
    lb = reshape(lb', [], 1);
    ub = reshape(ub', [], 1);

    % Reformulate the problem for ResQPASS
    % Try cholesky, else do svd decomposion
    [R,p] = chol(Q);
    if p ~= 0
        [V,D] = eig(Q);
        R = sqrt(max(D,0)) * V';
    end
    
    if preconditioning          % Determine preconditioner using iLU
        setup.type = 'ilutp';
        setup.droptol = 1e-3;
        [L,U] = ilu(sparse(Q), setup);
        M1 = @(x) U\(L\x);
    else
        M1 = [];
    end

    [sol, history] = ResQPASS_altered(R, zeros(Ndof, 1), lb, ub, ...
        "M1", M1);                  % Solve the contact problem
    
    % Update the points with the optimal differences
    sol_nodes = nodes + reshape(sol, 2, nV)';

    if preconditioning
        filename = filename + "_precond";
    end

    save("Create_2D_meshes/Results/" + filename + ".mat", "nodes", ...
        "triangles", "sol_nodes","sol", "history", "Q", "R", ...
        "lb_x", "ub_x", "lb_y", "ub_y");
end