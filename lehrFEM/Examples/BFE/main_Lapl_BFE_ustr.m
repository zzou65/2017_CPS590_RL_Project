% Run script for bilinear finite element solver.

% Copyright 2005-2005 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  theta = linspace(0,2*pi,7);
  theta = theta(1:6)';
  
  % Initialize constants
  
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);                                                       % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;                                                       % Dirichlet boundary data
  DHANDLE = inline(['dist_diff(dist_union(dist_tri(x,[0 1],[-sqrt(3)/2 -1/2],[sqrt(3)/2 -1/2]),' ...  % Signed distance function
                    'dist_tri(x,[0 -1],[sqrt(3)/2 1/2],[-sqrt(3)/2 1/2])),' ...
                    'dist_circ(x,[0 0],.5))'],'x');                                                                                                                                                          % Element size function
  FIXEDPOS = [cos(theta)/sqrt(3) sin(theta)/sqrt(3); ...                                              % Fixed positions of the mesh                                                                                                                          
              cos(theta+pi/6)    sin(theta+pi/6)];                                        
  
  % Initialize mesh
  
  Mesh = init_Mesh([-1 -1; 1 1],0.05,DHANDLE,@h_uniform,FIXEDPOS,0);
  Mesh = morph(Mesh);     
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
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
  