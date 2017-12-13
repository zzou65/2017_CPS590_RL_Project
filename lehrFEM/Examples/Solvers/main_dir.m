% Run script for direct solver.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  NREFS = 9;                                     % Number of red refinements
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data 
  
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
  L = assemLoad_LFE(Mesh,P3O3(),F_HANDLE);
    
  % Incorporate Dirichlet boundary conditions
  
  [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
  L = L - A*U;
  
  % Start direct solver
 
  t = cputime;
  U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
  fprintf('Runtime of direct solver [s]  :  %.2e\n',cputime-t);
  
  % Clear memory
  
  clear all;
  