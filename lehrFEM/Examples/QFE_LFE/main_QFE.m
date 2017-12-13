% Run script for piecewise quadratic finite elements 

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS =4;               % Number of red refinement steps
  F_HANDLE = @f_LShap;     % Right hand side source term
  GD_HANDLE = @g_D_LShap;  % Dirichlet boundary data
  GN_HANDLE = @g_N_LShap;  % Neumann boundary data
  
  NREFS =1;
      
  % Initialize mesh
  
  clear Mesh
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
 
  A = assemMat_QFE(Mesh,@STIMA_Lapl_QFE);
  L = assemLoad_QFE(Mesh,P7O6(),F_HANDLE);
   
  % Incorporate Neumann boundary data
  
  L = assemNeu_QFE(Mesh,-1:-1:-4,L,gauleg(0,1,4),GN_HANDLE);
  
  % Incorporate Dirichlet boundary data
 
  [U,FreeDofs] = assemDir_QFE(Mesh,-7,GD_HANDLE);
  L = L - A*U;
  
  % Solve the linear system
 
  U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
   
  % Plot out solution
  
  plot_QFE(U,Mesh);
  colorbar;
  plotLine_QFE(U,Mesh,[0,0.6],[1,0.6]);
  
  % Clear memory
  
  %clear all;
  