% Runs script for convergence rates of hp-FEM. This is for testing
% purposes.

% Copyright 2009 Christoph Wiesmeyr
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 15;                           % Number of mesh refinements  
  
%   F = @(x,varargin)zeros(size(x,1),1);  % Right hand side source term
%   GD = @gD_LShap;                     % Dirichlet boundary data
%   GRAD_UEX = @grad_uex_LShap;          % Gradient of exact solution
%   UEX = @uex;
  
  deg=5;
  F = @(x,varargin)-deg*(deg-1)*(x(:,1).^(deg-2).*x(:,2).^deg+x(:,2).^(deg-2).*x(:,1).^deg);  % Right hand side source term                       
  UEX = @(x,varargin)x(:,1).^deg.*x(:,2).^deg;
  GD = UEX;
  GRAD_UEX = @(x,varargin)deg*[x(:,2).^deg.*x(:,1).^(deg-1) x(:,2).^(deg-1).*x(:,1).^deg];
  
  
%   UEX = @(x,varargin)x(:,1).*x(:,2);
%   GD = UEX;
%   GRAD_UEX = @(x,varargin)[x(:,2) x(:,1)];
  
  % Initialize mesh

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;  
  
  CNodes = 3;
  
  % Prepare mesh for longest edge bisection 
  
  nDofs = zeros(1,NREFS);
  H1S_error = zeros(1,NREFS);
  H1S_errorEX = zeros(1,NREFS);
  Mesh = init_LEB(Mesh);
  for i = 1:NREFS
  
    Mesh = refine_hp(Mesh,CNodes);
    
    % Generate mesh data structure for hp-FEM
  
    Mesh_hp.Coordinates = Mesh.Coordinates;
    Mesh_hp.Elements = Mesh.Elements;
    Mesh_hp.ElemFlag = zeros(size(Mesh_hp.Elements,1),1);
    Mesh_hp = add_Edges(Mesh_hp);
    if(isfield(Mesh_hp,'BdFlags'))
      Mesh_hp = rmfield(Mesh_hp,'BdFlags');
    end
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
    
    % Compute discretization errors
    
    nDofs(i) = size(U,1);
    H1S_error(i) = H1SErr_hp(Mesh_hp,U,Elem2Dof,QuadRule_2D,Shap_2D,GRAD_UEX);
    
    % compute L^2 projection
    uint = comp_interp(Mesh_hp,CDofs,EDofs,QuadRule_2D,Shap_2D,UEX);
    H1S_errorEX(i) = H1SErr_hp(Mesh_hp,uint,Elem2Dof,QuadRule_2D,Shap_2D,GRAD_UEX);
    
    H1S_error(i)
    H1S_errorEX(i)
    
  end
  
  Errorplot(U,Mesh_hp,Elem2Dof,pmax,UEX); colorbar;
    
  % Generate figure
  
  fig = figure('Name','Convergence rates for hp-FEM');
  plot(sqrt(nDofs),H1S_error,'r-o');
  title('{\bf Discretization errors for hp-FEM}')
  xlabel('{\bf # Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  set(gca,'XScale','lin','YScale','log');
 
  % Clear memory
  
  %clear all;
  