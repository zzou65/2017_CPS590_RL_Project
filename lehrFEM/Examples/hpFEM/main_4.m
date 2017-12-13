% Runs script for convergence rates of hp-FEM.

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  PMAX = 10;                            % Maximum polynomial degree
  F = @(x,varargin)zeros(size(x,1),1);  % Right hand side source term
  GD = @gD_LShap;                       % Dirichlet boundary data
  GRAD_UEX = @grad_uex_LShap;           % Gradient of exact solution
  
  % Initialize mesh

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;  
  Mesh = add_Edge2Elem(Mesh);
  
  
  nDofs = zeros(1,PMAX);
  H1S_error = zeros(1,PMAX);
  for p = 1:PMAX
    
    % Assign polynomial degrees and build dof maps
  
    EDofs = (p-1)*ones(size(Mesh.Edges,1),1);
    if(p > 2)
      CDofs = (p-1)*(p-2)/2*ones(size(Mesh.Elements,1),1);
    else
      CDofs = zeros(size(Mesh.Elements,1),1);
    end
    Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs);
  
    % Build shape functions and quadrature rules

    QuadRule_1D = gauleg(0,1,2*p);
    Shap_1D = shap_hp([QuadRule_1D.x zeros(size(QuadRule_1D.x))],p);
    QuadRule_2D = Duffy(TProd(QuadRule_1D));
    Shap_2D = shap_hp(QuadRule_2D.x,p);
  
    % Assemble global load vector and mass matrix
  
    A = assemMat_hp(Mesh,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D);
    L = assemLoad_hp(Mesh,Elem2Dof,QuadRule_2D,Shap_2D,F);
  
    % Incoporate Dirichlet boundary conditions
  
    [U,FreeDofs] = assemDir_hp(Mesh,Elem2Dof,-1,QuadRule_1D,Shap_1D,GD);
    L = L - A*U;
  
    % Solve the linear system
  
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        
    % Compute discretization errors
    
    nDofs(p) = size(U,1);
    H1S_error(p) = H1SErr_hp(Mesh,U,Elem2Dof,QuadRule_2D,Shap_2D,GRAD_UEX);
    
  end
    
  % Generate figure
  
  fig = figure('Name','Convergence rates for hp-FEM');
  plot(nDofs,H1S_error,'r-o');
  title('{\bf Discretization errors for hp-FEM}')
  xlabel('{\bf # Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  set(gca,'XScale','log','YScale','log');
 
  % Clear memory
  
  clear all;
  