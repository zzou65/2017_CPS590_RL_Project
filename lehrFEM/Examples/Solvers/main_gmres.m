% Run script for GMRES solver.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  NREFS = 7;                                     % Number of red refinements
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data 
  TOL = 1e-12;                                   % Stopping criterion
  MAXIT = 100;                                   % Maximum number of iterations
  RESTART = 10;                                  % Restart parameter 
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
      
  % Generate multigrid data structure
  
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);   
  end
    
  % Compute mass matrix and load vector
  
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  L = assemLoad_LFE(Mesh,P1O2(),F_HANDLE);
    
  % Incorporate Dirichlet boundary conditions
  
  [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
  L = L - A*U;
  A = A(FreeDofs,FreeDofs);
  L = L(FreeDofs);
  
  % Start mulitgrid solver
 
  t = cputime;
  U0 = zeros(size(L));
  [U(FreeDofs),flag,relres,iter,resvec] = gmres_solve(U0,A,L, ...
                                                      RESTART,TOL,MAXIT);
  fprintf('Runtime of GMRES solver [s]  :  %f\n',cputime-t);
    
  fig = figure('Name','GMRES solver');
  plot(resvec,'rx');
  title('{\bf GMRES solver}');
  xlabel('{\bf Iteration number}');
  ylabel('{\bf Relative residual [log]}');
  set(gca,'YScale','log');
  
  % Clear memory
  
  clear all;
  