% Convergence rates for piecewise quadratic and linear finite elements for
% the Poisson equation with Dirichlet boundary conditions on the unit
% ball. This script generates the following .eps files:
%  B2_meshwidth.eps,  B2_dofs.eps.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  DHANDLE = @dist_circ;   % Signed distance function
  C = [0 0];              % Center of the circle
  R = 1;                  % Radius of the circle
  INIT_NREFS = 2;         % Number of initial mesh refinemenets
  NREFS = 4;              % Number of red refinement steps
  F_HANDLE = @f_Ball;     % Right hand side source term
  GD_HANDLE = @g_D_Ball;  % Dirichlet boundary data
  MAXIT = 41;             % Maximum number of CG iterations
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Do initial mesh refinements
  
  for i = 1:INIT_NREFS
    Mesh = refine_REG(Mesh,DHANDLE,C,R);
  end
  
  % Compute discretization error on a series of meshes
    
  rho_LFE = zeros(NREFS,MAXIT);
  nDofs_LFE = zeros(1,NREFS);
  rho_QFE = zeros(NREFS,MAXIT);
  nDofs_QFE = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh,DHANDLE,C,R);    
    
    % Assemble Stiffness matrix, load vector and incorporate BC
  
    A_QFE = assemMat_QFE(Mesh,@STIMA_Lapl_QFE);
    A_LFE = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    
    % Discard Dirichlet nodes from matrices

    [U_LFE,FreeDofs_LFE] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    A_LFE = A_LFE(FreeDofs_LFE,FreeDofs_LFE);
    U_LFE = U_LFE(FreeDofs_LFE);
    
    [U_QFE,FreeDofs_QFE] = assemDir_QFE(Mesh,-1,GD_HANDLE);
    A_QFE = A_QFE(FreeDofs_QFE,FreeDofs_QFE);
    U_QFE = U_QFE(FreeDofs_QFE);
    
    % Solve the linear system
  
    nDofs_LFE(i) = size(U_LFE,1);
    x0 = zeros(size(U_LFE));
    b = ones(size(U_LFE));
    [x,flag,relres,iter,resvec] = cg_solve(x0,A_LFE,b,0,MAXIT);
    rho_LFE(i,:) = resvec;
    
    nDofs_QFE(i) = size(U_QFE,1);
    x0 = zeros(size(U_QFE));
    b = ones(size(U_QFE));
    [x,flag,relres,iter,resvec] = cg_solve(x0,A_QFE,b,0,MAXIT);
    rho_QFE(i,:) = resvec;
    
  end
  
  rate_LFE = (rho_LFE(:,MAXIT) ./ rho_LFE(:,1)).^(1/(MAXIT-1))
  rate_QFE = (rho_QFE(:,MAXIT) ./ rho_QFE(:,1)).^(1/(MAXIT-1))
  
  % Generate figures

  fig = figure();
  cg_it = 0:(MAXIT-1);
  plot(cg_it,rho_LFE(1,:),'r-', ...
       cg_it,rho_LFE(2,:),'b-', ...
       cg_it,rho_LFE(3,:),'g-', ...
       cg_it,rho_LFE(4,:),'m-');
  xlabel('{\bf CG iterations}');
  ylabel('{\bf relative Euclidean residual}');
  legend(sprintf('N = %d',nDofs_LFE(1)), ...
         sprintf('N = %d',nDofs_LFE(2)), ...
         sprintf('N = %d',nDofs_LFE(3)), ...
         sprintf('N = %d',nDofs_LFE(4)), ...
         'Location','SouthWest');
  set(gca,'YScale','log');
  
  print('-depsc','CG_rate_Ball.eps');
  close(fig);
  !gv CG_rate_Ball.eps &
  
