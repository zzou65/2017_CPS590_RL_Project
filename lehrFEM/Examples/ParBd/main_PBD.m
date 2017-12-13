% Run script for finite element solver with parabolic boundary
% approximation. 
 
%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 5;             % Number of inital mesh refinements
  DHANDLE = @dist_circ;  % Signed distance function
  C = [0 0];             % Center of the circle
  R = 1;                 % Radius of the circle
  F_HANDLE = @f;         % Right hand side source term
  GD_HANDLE = @g_D;      % Dirichlet boundary data
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  for i = 1:NREFS
    Mesh = refine_REG(Mesh,DHANDLE,C,R);    
  end  
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_ParBd(Mesh,DHANDLE,C,R);
    
  % Assemble mass matrix and load vector
 
  A_PBD = assemMat_PBD(Mesh,@STIMA_Lapl_PBD,P7O6());
  L_PBD = assemLoad_PBD(Mesh,P7O6(),F_HANDLE);
    
  % Incorporate Dirichlet boundary conditions
  
  [U_PBD,FreeDofs] = assemDir_PBD(Mesh,-1,GD_HANDLE);
  L_PBD = L_PBD - A_PBD*U_PBD;
    
  % Solve the linear system
  
  U_PBD(FreeDofs) = A_PBD(FreeDofs,FreeDofs)\L_PBD(FreeDofs);
      
  % Generate figures
  
  plot_PBD(U_PBD,Mesh);
  colorbar;
  title('{\bf Quadratic finite elements with parabolic boundaries}');
  
  % Clear memory
  
  clear all;
  