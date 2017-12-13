% Run script for preconditioned GMRES solver.

%   Copyright 200-2006 Patrick Meury
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
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  CMesh = add_Edges(CMesh);
  Loc = get_BdEdges(CMesh);
  CMesh.BdFlags = zeros(size(CMesh.Edges,1),1);
  CMesh.BdFlags(Loc) = -1;
  CMesh.ElemFlag = zeros(size(CMesh.Elements,1),1);
  
  % Compute non-Dirichlet vertices
  
  Loc = get_BdEdges(CMesh);
  DDofs = unique([CMesh.Edges(Loc,1); CMesh.Edges(Loc,2)]);
  CFDofs = setdiff(1:size(CMesh.Coordinates,1),DDofs);
  
  % Generate multigrid data structure
  
  for i = 1:NREFS
        
    % Compute stiffness matrix and load vector
    
    A = assemMat_LFE(CMesh,@STIMA_Lapl_LFE);
    ML_Data.D{i} = diag(A(CFDofs,CFDofs));
    
    % Refine the mesh and compute prolongation matrix
    
    FMesh = refine_REG(CMesh);   
    P = get_PMat_LFE(CMesh,FMesh);
    Loc = get_BdEdges(FMesh);
    DDofs = unique([FMesh.Edges(Loc,1); FMesh.Edges(Loc,2)]);
    FFDofs = setdiff(1:size(FMesh.Coordinates,1),DDofs);    
    ML_Data.P{i} = P(FFDofs,CFDofs);
    
    % Update coarse mesh
    
    CMesh = FMesh;
    CFDofs = FFDofs;
    
  end
  
  % Compute load vector on finest mesh
  
  A = assemMat_LFE(CMesh,@STIMA_Lapl_LFE);
  L = assemLoad_LFE(CMesh,P1O2(),F_HANDLE);
  
  % Incoporate Dirichlet boundary conditions
  
  [U,CFDofs] = assemDir_LFE(CMesh,-1,GD_HANDLE);
  L = L - A*U;
  A = A(CFDofs,CFDofs);
  ML_Data.D{NREFS+1} = diag(A);
  L = L(CFDofs);
  
  % Run preconditioned CG solver
 
  t = cputime;
  U0 = L;
  [U(CFDofs),flag,relres,iter,resvec] = pgmres_solve(U0,A,L, ...
                                                     RESTART,TOL,MAXIT, ...
                                                     @bpx_prec,ML_Data);
  fprintf('Runtime of preconditioned GMRES solver [s]  :  %f\n',cputime-t);
   
  fig = figure('Name','Preconditioned GMRES solver');
  plot(resvec,'rx');
  title('{\bf Preconditioned GMRES solver (BPX preconditioner)}');
  xlabel('{\bf Iteration number}');
  ylabel('{\bf Relative residual [log]}');
  set(gca,'YScale','log');
  
  % Clear memory
  
  clear all;
  