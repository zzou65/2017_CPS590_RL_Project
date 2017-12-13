% Generates a quadrilateral structured mesh of L-shape with Dirichlet and 
% Neumann boundary conditions

% Copyright 2005-2005 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
 
  % Initialize constants
   
  NREFS = 4;                  % Number of red refinement steps
  F_HANDLE = @f_LShap;        % Right hand side source term
  GD_HANDLE = @g_D_LShap;     % Dirichlet boundary data
  GN_HANDLE = @g_N_LShap;     % Neumann boundary data

  % Initialize mesh
   
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = [-1 -2 -7 -7 -3 -4 -4 -2];
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end    
  Mesh = add_Edge2Elem(Mesh);
    
  % Assemble stiffness matrix and load vector
 
  QuadRule = TProd(gauleg(0,1,2));
  A = assemMat_BFE(Mesh,@STIMA_Lapl_BFE,QuadRule);
  L = assemLoad_BFE(Mesh,QuadRule,F_HANDLE);
  
  % Incorporate Neumann boundary data
  
  L = assemNeu_BFE(Mesh,-1:-1:-2,L,gauleg(0,1,2),GN_HANDLE);
   
  % Incorporate Dirichlet boundary data
 
  [U,FreeDofs] = assemDir_BFE(Mesh,[-3 -4 -7],GD_HANDLE);
  L = L - A*U;
  
  % Solve the linear system
 
  U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
  % Plot out solution
    
  plot_BFE(U,Mesh);
  colorbar;
  
  % Clear memory
  
  clear all;
  