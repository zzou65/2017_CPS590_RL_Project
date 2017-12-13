% Run script for CG solver with Bramble-Pasciak-Xu preconditioner.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 6;                                     % Number of red refinements
  U_Handle = @(x,varargin)ones(size(x,1),1);
  F_Handle = @(x,varargin)pi^2*[sin(pi*x(:,2)) sin(pi*x(:,1))]+[sin(pi*x(:,2)) sin(pi*x(:,1))];
  GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
  TOL = 1e-12;                                   % Stopping criterion
  MAXIT = 2000;                                  % Maximum number of iterations
    
  % Initialize mesh
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  CMesh = add_Edges(CMesh);
  Loc = get_BdEdges(CMesh);
  CMesh.BdFlags = zeros(size(CMesh.Edges,1),1);
  CMesh.BdFlags(Loc) = -1;
  CMesh.ElemFlag = zeros(size(CMesh.Elements,1),1);
  
  % Compute non-Dirichlet vertices
  
  Loc = get_BdEdges(CMesh);
  DNodes = unique([CMesh.Edges(Loc,1); CMesh.Edges(Loc,2)]);
  CFNodes = setdiff(1:size(CMesh.Coordinates,1),DNodes);
  nCCoordinates = size(CMesh.Coordinates,1);
  CFDofs = [CFNodes;CFNodes+nCCoordinates];

  % Generate multigrid data structure
  
  for i = 1:NREFS
    
    
    % Compute stiffness matrix and load vector
    
    [IC,JC,C] = assemMat_LFE2(CMesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
    [IM,JM,M] = assemMat_LFE2(CMesh,@MASS_LFE2);
    A = sparse([IC;IM],[JC;JM],[C;M]);
    ML_Data.D{i} = diag(A(CFDofs,CFDofs));
    
    % Refine the mesh and compute prolongation matrix
    
    FMesh = refine_REG(CMesh);   
    nFCoordinates = size(FMesh.Coordinates,1);
    P = get_PMat_LFE2(CMesh,FMesh);
    Loc = get_BdEdges(FMesh);
    DNodes = unique([FMesh.Edges(Loc,1); FMesh.Edges(Loc,2)]);
    FFNodes = setdiff(1:size(FMesh.Coordinates,1),DNodes);    
    FFDofs = [FFNodes;FFNodes+nFCoordinates];
    ML_Data.P{i} = P(FFDofs,CFDofs);
    
    % Update coarse mesh
    
    CMesh = FMesh;
    CFDofs = FFDofs;
    
  end
  
  % Compute load vector on finest mesh
  
  [IC,JC,C] = assemMat_LFE2(CMesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
  [IM,JM,M] = assemMat_LFE2(CMesh,@MASS_LFE2);
  A = sparse([IC;IM],[JC;JM],[C;M]);
  L = assemLoad_LFE2(CMesh,P7O6(),F_Handle);
  
  % Incoporate Dirichlet boundary conditions
  
  [U,CFDofs] = assemDir_LFE2(CMesh,-1,GD_Handle);
  L = L - A*U;
  A = A(CFDofs,CFDofs);
  ML_Data.D{NREFS+1} = diag(A);
  L = L(CFDofs);
  
  % Run preconditioned CG solver
 
  t = cputime;
  U0 = L;
  [U(CFDofs),flag,relres,iter,resvec] = pcg_solve(U0,A,L,TOL,MAXIT, ...
                                                  @bpx_prec,ML_Data);
  fprintf('Runtime of preconditioned CG solver [s]  :  %f\n',cputime-t);
  
  fig = figure('Name','Preconditioned CG solver');
  plot(resvec,'rx');
  title('{\bf Preconditioned CG solver (BPX preconditioner)}');
  xlabel('{\bf Iteration number}');
  ylabel('{\bf Relative residual [log]}');
  set(gca,'YScale','log');
  
  % Clear memory
  
  clear all;
  
  