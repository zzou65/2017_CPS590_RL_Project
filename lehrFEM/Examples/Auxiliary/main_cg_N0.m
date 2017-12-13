% Run script for the auxiliary preconditioner to slove edge element using
% N0 solution in Nodal element space (with complete Dirichlet boundary conditions)

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
    clear Mesh;
    
    % Initialize constant

    NREFS = 5;
    Tal = .5;
    TOL = 1e-12;                                   % Stopping criterion
    MAXIT = 2000;                                  % Maximum number of iterations
    NTimes = 2;
    U_Handle = @(x,varargin) ones(size(x,1),1);
    F_Handle = @(x,varargin) (pi^2+Tal)*[sin(pi*x(:,2)) sin(pi*x(:,1))];
    GD_Handle = @(x,varargin) [sin(pi*x(:,2)) sin(pi*x(:,1))];
%     F_Handle = @(x,varargin)Tal*[-x(:,2) x(:,1)];
%     GD_Handle = @(x,varargin)[-x(:,2) x(:,1)];
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    for i=1:NREFS
        Mesh = refine_REG(Mesh);
    end
    
    % Assemble stiffness matrix, Curl-curl matrix, MASS matrix and load vector

    % A = stiffness matrix for Laplacian (pure Neumann problem)
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    
    % Assemble stiffness matrix for curl*\mu^{-1}*curl-operator
    % Function \mu can be passed through U_HANDLE
    [ICe,JCe,Ce] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,U_Handle,P7O6());
    % Assemble mass matrix for edge elements
    % TODO: General coefficient function
    [IMe,JMe,Me] = assemMat_W1F(Mesh,@MASS_W1F,U_Handle,P7O6());
    % Build complete edge element Galerkin matrix Se
    Se = sparse([ICe;IMe],[JCe;JMe],[Ce;Tal*Me]);
    % Right hand side vector for edge elements
    Le = assemLoad_W1F(Mesh,P7O6(),F_Handle);
    
    % Assemble matrix for vector Laplacian with respect to 
    % continuous linear Lagrangian finite elements
    [ICn,JCn,Cn] = assemMat_LFE2(Mesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
    [IDn,JDn,Dn] = assemMat_LFE2(Mesh,@STIMA_Div_LFE2,U_Handle,P7O6());
    [IMn,JMn,Mn] = assemMat_LFE2(Mesh,@MASS_LFE2);
    Sn = sparse([ICn;IDn;IMn],[JCn;JDn;JMn],[Cn;Dn;Tal*Mn]);
    
    % Prolongation matrix (nodal representation -> edge element representation)
    P = convert_LFE2_W1F(Mesh);
    
    % Gradient matrix (edge-vertex incidence matrix)
    GE = Mat_G2E(Mesh);
    
    % Incorporate Dirichlet boundary data
    
    Loc = get_BdEdges(Mesh);
    BdNodes = unique(Mesh.Edges(Loc,:));
    FreeDofs_LFE = setdiff(1:size(Mesh.Coordinates),BdNodes);
    
    [Ue,FreeDofs_e] = assemDir_W1F(Mesh,-1,GD_Handle,gauleg(0,1,1));
    Le = Le - Se*Ue;
    [Un,FreeDofs_n] = assemDir_LFE2(Mesh,-1,GD_Handle);
    
    Se = Se(FreeDofs_e,FreeDofs_e);
    Le = Le(FreeDofs_e);
    Sn = Sn(FreeDofs_n,FreeDofs_n);
    P = P(FreeDofs_e,FreeDofs_n);
    A = A(FreeDofs_LFE,FreeDofs_LFE);
    inv_A = inv(A);
    GE = GE(FreeDofs_e,FreeDofs_LFE);
    T_GE = transpose(GE);
    
    % Run pcg solver
    
    t = cputime;
    U0 = Le;
    [Ue(FreeDofs_e),flag,relres,iter,resvec] = pcg_solve(U0,Se,Le,TOL,MAXIT,@aux_prec,...
        Se,P,NTimes,inv_A,GE,T_GE,@solve_N0,Sn);
    fprintf('Runtime of preconditioned CG solver [s]  :  %f\n',cputime-t);

    % Plot figures
    
    fig = figure('Name','Preconditioned CG solver');
    plot(resvec,'rx');
    title('{\bf Preconditioned CG solver (N0 solution in nodal FE)}');
    xlabel('{\bf Iteration number}');
    ylabel('{\bf Relative residual [log]}');
    set(gca,'YScale','log');
    grid on
    
    % Output .eps files
    
    print('-depsc', 'CG_prec_N0.eps');
    !gv CG_prec_N0.eps &
    