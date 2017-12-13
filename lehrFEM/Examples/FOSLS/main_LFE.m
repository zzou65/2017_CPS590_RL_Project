% Run script for piecewise linear finite element solver.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
  NREFS = 5;               % Number of red refinement steps
  F_HANDLE = @f_LShap;     % Right hand side source term
  GD_HANDLE = @g_D_LShap;  % Dirichlet boundary data
  GN_HANDLE = @g_N_LShap;  % Neumann boundary data
 
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = [-1 -2 -7 -7 -3 -4];
  for i = 1:NREFS
   Mesh = refine_REG(Mesh);     
  end
  Mesh = add_Edge2Elem(Mesh);
  
  % Assemble stiffness matrix and load vector
 
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE,0,1,2);
   
  % Incorporate Neumann boundary data
  
  L = assemNeu_LFE(Mesh,-1:-1:-4,L,gauleg(0,1,4),GN_HANDLE);
  
  % Incorporate Dirichlet boundary data
 
  [U,FreeDofs] = assemDir_LFE(Mesh,-7,GD_HANDLE);
  L = L - A*U;
  
  % Solve the linear system
 
  U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
  % Plot out solution
    
  plot_LFE(U,Mesh);
  colorbar;
  
  % Clear memory
  
  clear all;
  