% Runs script for hp-FEM. convergence bilinear value

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  clear Mesh;
  
  NREFS = 10;                                                         % Number of mesh refinements  
  F = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD = @(x,varargin)zeros(size(x,1),1);                       % Dirichlet boundary data
  v_Handle = @(x,varargin)ones(size(x,1),2);               % Dirichlet boundary data
  
  U=@(x,varargin) x(:,1).^2+x(:,2).^2;
  CU=@(x,varargin) 2*x(:,1)+2*x(:,2);
  V=@(x,varargin) x(:,1)+x(:,2);

  % Initialize mesh

  Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;  
  
  CNodes = transpose(1:size(Mesh.Coordinates,1));
  
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
  
  [EDofs,CDofs] = assign_pdeg(Mesh_hp,CNodes,NREFS);
%   EDofs=zeros(size(Mesh_hp.Edges,1),1);
%   CDofs=zeros(size(Mesh_hp.Elements,1),1);
  Elem2Dof = build_DofMaps(Mesh_hp,EDofs,CDofs);
  pmax = max(CDofs)+1;
  
  % Build shape functions and quadrature rules

  QuadRule_1D = gauleg(0,1,2*pmax);
  Shap_1D = shap_hp([QuadRule_1D.x zeros(size(QuadRule_1D.x))],pmax);
  QuadRule_2D = Duffy(TProd(QuadRule_1D));
  Shap_2D = shap_hp(QuadRule_2D.x,pmax);
  velocity=v_Handle(QuadRule_2D.x);
  
  % Assemble global load vector and mass matrix
  
  %A = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D);
  M = assemMat_hp(Mesh_hp,Elem2Dof,@MASS_hp,QuadRule_2D,Shap_2D);
  B = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Conv_hp,QuadRule_2D,Shap_2D,velocity);
  LU = assemLoad_hp(Mesh_hp,Elem2Dof,QuadRule_2D,Shap_2D,U);
  LV = assemLoad_hp(Mesh_hp,Elem2Dof,QuadRule_2D,Shap_2D,V);
  
  u=M\LU;
  v=M\LV;
  
  norm(u'*B*v-innerproduct(CU,V,-1,1,-1,1))
  % Incoporate Dirichlet boundary conditions
  
  %[U,FreeDofs] = assemDir_hp(Mesh_hp,Elem2Dof,-1,QuadRule_1D,Shap_1D,GD);
  %L = L - A*U;
  
  % Solve the linear system
  
  %U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
 
  % Plot hp-FEM solution
 
  plot_hp(u,Mesh_hp,Elem2Dof,pmax);
  colorbar;
  plot_hp(B*u,Mesh_hp,Elem2Dof,pmax);
  colorbar;
    
  fig = figure('Name','Polynomial degrees');  
  patch('Faces',Mesh_hp.Elements, ...
        'Vertices',Mesh_hp.Coordinates, ...
        'FaceVertexCData',CDofs, ...
        'EdgeColor','k', ...
        'FaceColor','flat');
  set(gca,'CLim',[1 NREFS],'DataAspectRatio',[1 1 1]);
  colormap(jet);
  alpha(.9);
  colorbar;
  set(gcf,'renderer','openGL');
  % Clear memory
  
  clear all;