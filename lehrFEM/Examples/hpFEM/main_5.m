% Runs script for hp-FEM.

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 18;                                               % Number of mesh refinements  
  F = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD = @(x,varargin)zeros(size(x,1),1);                     % Dirichlet boundary data
  
  % Initialize mesh

  Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;  
  
  % CNodes = transpose(1:size(Mesh.Coordinates,1));
  CNodes = 1;
  
  % Prepare mesh for longest edge bisection 
  
  Mesh = init_LEB(Mesh);
  for i = 1:NREFS
    Mesh = refine_hp(Mesh,CNodes);
  end 
  
  % Generate mesh data structure for hp-FEM
  
  Mesh_hp.Coordinates = Mesh.Coordinates;
  Mesh_hp.Elements = Mesh.Elements;
  Mesh_hp.ElemFlag = zeros(size(Mesh_hp.Elements,1),1);
  Mesh_hp = add_Edges(Mesh_hp);
  Loc = get_BdEdges(Mesh_hp);
  Mesh_hp.BdFlags = zeros(size(Mesh_hp.Edges,1),1);
  Mesh_hp.BdFlags(Loc) = -1;
  Mesh_hp = add_Edge2Elem(Mesh_hp);
  
  % Assign polynomial degrees and build dof maps
  
  [EDofs,CDofs,ElemDeg] = assign_pdeg(Mesh_hp,CNodes,NREFS);
  Elem2Dof = build_DofMaps(Mesh_hp,EDofs,CDofs);
  pmax = max(ElemDeg);
  
  % Build shape functions and quadrature rules

  QuadRule_1D = gauleg(0,1,2*pmax);
  Shap_1D = shap_hp([QuadRule_1D.x zeros(size(QuadRule_1D.x))],pmax);
  QuadRule_2D = Duffy(TProd(QuadRule_1D));
  Shap_2D = shap_hp(QuadRule_2D.x,pmax);
  
  % Assemble global load vector and mass matrix
  
  A = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D);
  L = assemLoad_hp(Mesh_hp,Elem2Dof,QuadRule_2D,Shap_2D,F);
  
  % Incoporate Dirichlet boundary conditions
  
  [U,FreeDofs] = assemDir_hp(Mesh_hp,Elem2Dof,-1,QuadRule_1D,Shap_1D,GD);
  L = L - A*U;
  
  % Solve the linear system
  
  U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
 
  % Plot hp-FEM solution
 
  plot_hp(U,Mesh_hp,Elem2Dof,pmax);
  colorbar;
    
  % Clear memory
  
  clear all;
  