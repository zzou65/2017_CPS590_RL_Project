% Run script for bilinear finite element solver.

% Copyright 2005-2005 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
 
  % Initialize constants
  
  NREFS = 3;               % Number of unifrom red refinements
  F_HANDLE = @f_LShap;     % Right hand-side source term 
  GD_HANDLE = @g_D_LShap;  % Dirichlet boundary data 
  
  % Initialize mesh
   
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);   
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end    
        
  % Assemble stiffness matrix and load vector
 
  QuadRule = TProd(gauleg(0,1,2));
  A = assemMat_BFE(Mesh,@STIMA_Lapl_BFE,QuadRule);
  L = assemLoad_BFE(Mesh,QuadRule,F_HANDLE);
   
  % Incorporate Dirichlet boundary data
 
  [U,FreeDofs] = assemDir_BFE(Mesh,-1,GD_HANDLE);
  L = L - A*U;
  
  % Solve the linear system
 
  U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
  % Plot out solution
    
  plot_BFE(U,Mesh);
  colorbar;
  
  % Clear memory
  
  clear all;
  