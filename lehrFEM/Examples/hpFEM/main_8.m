% run script for convergence rates of p-fem

% Copyright 2009 Christoph Wiesmeyr
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  clear;


  % Initialize constants
  
  NREFS = 15;                                                     % Number of mesh refinements  
  
  F = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));         % Right hand side source term
  GD = @(x,varargin)zeros(size(x,1),1);                           % Dirichlet boundary data
  UEX = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));
  GRAD_UEX = @(x,varargin)pi*[cos(pi*x(:,1)).*sin(pi*x(:,2)) ...  % Gradient of exact solution
                              sin(pi*x(:,1)).*cos(pi*x(:,2))]; 
                          
 
  % Initialize mesh

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;  
  Mesh = add_Edge2Elem(Mesh);
  
  H1S_error = zeros(1,NREFS);
  H1S_errorEX = zeros(1,NREFS);
  
  for i = 1:NREFS 
      
    pmax=i; % do p-refinement
    
    EDofs = (pmax-1)*ones(size(Mesh.Edges,1),1);
    if(pmax > 2)
      CDofs = (pmax-1)*(pmax-2)/2*ones(size(Mesh.Elements,1),1);
    else
      CDofs = zeros(size(Mesh.Elements,1),1);
    end

    Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs);  
   
      
    % Build shape functions and quadrature rules
    QuadRule_1D = gauleg(0,1,2*pmax);
    Shap_1D = shap_hp([QuadRule_1D.x zeros(size(QuadRule_1D.x))],pmax);
    QuadRule_2D = Duffy(TProd(QuadRule_1D));
    Shap_2D = shap_hp(QuadRule_2D.x,pmax);
  
    % Assemble global load vector and mass matrix
  
    A = assemMat_hp(Mesh,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D);
    L = assemLoad_hp(Mesh,Elem2Dof,QuadRule_2D,Shap_2D,F);
  
    % Incoporate Dirichlet boundary conditions
  
    [U,FreeDofs] = assemDir_hp(Mesh,Elem2Dof,-1,QuadRule_1D,Shap_1D,GD);
    L = L - A*U;
  
    % Solve the linear system
  
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        
    % Compute discretization errors
    
    nDofs(i) = size(U,1);
    H1S_error(i) = H1SErr_hp(Mesh,U,Elem2Dof,QuadRule_2D,Shap_2D,GRAD_UEX);
        
  end
  
  % Generate figure
  
  fig = figure('Name','Convergence rates for p-FEM');
  plot(sqrt(nDofs),H1S_error,'r-o');
  title('{\bf Discretization errors for p-FEM}')
  xlabel('{\bf # Dofs [sqrt]}');
  ylabel('{\bf Discretization error [log]}');
  set(gca,'XScale','lin','YScale','log');
 
  % Clear memory
  
  clear all;
  