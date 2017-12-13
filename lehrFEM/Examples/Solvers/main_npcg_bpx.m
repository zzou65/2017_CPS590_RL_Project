% Run script for nested CG solver with Bramble-Pasciak-Xu preconditioner.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  IREFS = 2;                                     % Initial red refinements (2 or greater)
  NREFS = 7;                                     % Number of red refinements
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-12;                                   % Stopping criterion
  MAXIT = 2000;                                  % Maximum number of iterations
    
  % Initialize mesh
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  CMesh = add_Edges(CMesh);
  Loc = get_BdEdges(CMesh);
  CMesh.BdFlags = zeros(size(CMesh.Edges,1),1);
  CMesh.BdFlags(Loc) = -1;
  CMesh.ElemFlag = zeros(size(CMesh.Elements,1),1);
  
  for i = 1:IREFS
    CMesh = refine_REG(CMesh);
  end
  
  % Set up prescribed tolerance
  
  tol = TOL*4^(IREFS+NREFS-1);
  
  % Generate multigrid data structure
 
  A = assemMat_LFE(CMesh,@STIMA_Lapl_LFE);
  L = assemLoad_LFE(CMesh,P3O3(),F_HANDLE);
  [U,CFDofs] = assemDir_LFE(CMesh,-1,GD_HANDLE);
  L = L - A*U;
  U0 = zeros(size(CFDofs,2),1);
  total = 0;
  for i = 1:NREFS 
      
    A = A(CFDofs,CFDofs);
    ML_Data.D{i} = diag(A);
    L = L(CFDofs);
    
    % Compute solution on the current level
    
    t0 = cputime();
    [U0,flag,relres,iter,resvec] = pcg_solve(U0,A,L,TOL,MAXIT, ...
                                             @bpx_prec,ML_Data);
    t1 = cputime()-t0;
    if(flag == 0)
      error('PCG solver failed to converge.');  
    else
      fprintf('\n');
      fprintf('Level %2d\n',i);
      fprintf('  Runtime [s]         :  %.1e\n',t1);
      fprintf('  Relative tolerance  :  %.1e\n',tol);
      fprintf('  PCG iterations      :  %d\n',iter);
      fprintf('\n');
    end
    
    % Refine the mesh
    
    FMesh = refine_REG(CMesh);
    
    % Compute stiffness matrix and load vector
    
    A = assemMat_LFE(FMesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(FMesh,P1O2(),F_HANDLE);
    [U,FFDofs] = assemDir_LFE(FMesh,-1,GD_HANDLE);
    L = L-A*U;
    
    % Compute prolongation matrix
    
    P = get_PMat_LFE(CMesh,FMesh);
    P = P(FFDofs,CFDofs);
    ML_Data.P{i} = P;
    
    % Update coarse mesh
    
    CMesh = FMesh;
    CFDofs = FFDofs;
    
    % Prolongate initial value
    
    U0 = P*U0;
    
    % Update error tolerance
    
    tol = tol/4;
    
  end
  
  % Compute load vector on finest mesh
  
  A = assemMat_LFE(CMesh,@STIMA_Lapl_LFE);
  L = assemLoad_LFE(CMesh,P1O2(),F_HANDLE);
  
  % Incoporate Dirichlet boundary conditions
  
  [u,CFNodes] = assemDir_LFE(CMesh,-1,GD_HANDLE);
  L = L - A*U;
  A = A(CFDofs,CFDofs);
  ML_Data.D{NREFS+1} = diag(A); 
  L = L(CFDofs);
  
  % Run preconditioned CG solver
 
  t0 = cputime();
  [U(CFDofs),flag,relres,iter,resvec] = pcg_solve(U0,A,L,TOL,MAXIT, ...
                                                  @bpx_prec,ML_Data);
  t1 = cputime()-t0;
  total = total+t1;
  if(flag == 0)
    error('PCG solver failed to converge.');  
  else
    fprintf('\n');
    fprintf('Level %2d\n',NREFS+1);
    fprintf('  Runtime [s]         :  %.1e\n',t1);
    fprintf('  Relative tolerance  :  %.1e\n',tol);
    fprintf('  PCG iterations      :  %d\n',iter);
    fprintf('\n');
  end
    
  fprintf('Total runtime [s]  :  %.1e\n',total);
  
  % Generate figure
  
  fig = figure('Name','Preconditioned CG solver');
  plot(resvec,'rx');
  title('{\bf Preconditioned CG solver (BPX preconditioner)}');
  xlabel('{\bf Iteration number}');
  ylabel('{\bf Relative residual [log]}');
  set(gca,'YScale','log');
  
  % Clear memory
  
  clear all;
  