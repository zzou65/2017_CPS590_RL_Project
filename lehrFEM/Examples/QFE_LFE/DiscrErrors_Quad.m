% Convergence rates for piecewise quadratic finite elements for the Laplace
% equation with Dirichlet boundary conditions on the square for different
% quadrature rules.
 
%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 5;                                                      % Number of red refinement steps
  F_HANDLE = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));        % Dirichlet boundary data
  U_EX = @(x,varargin)pi*[sin(pi*x(:,2)).*cos(pi*x(:,1)) ...      % Exact solution for H1 semi-norm
                            sin(pi*x(:,1)).*cos(pi*x(:,2))];                                                                  
  SIGMA_HANDLE = @(x,varargin)ones(size(x,1),1);                  % Heat conductivity
  JIG = 1;                                                        % Jiggle parameter
    
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Compute discretization error on a series of meshes
    
  h = zeros(1,NREFS);
  N = zeros(1,NREFS);
  H1SErr_LO = zeros(1,NREFS);
  H1SErr_HI = zeros(1,NREFS); 
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh);    
    Mesh = add_Edge2Elem(Mesh);
    
    % Mesh preprocessing
    
    switch(JIG)
      case 1
        NewMesh = Mesh;      
      case 2
        Loc = get_BdEdges(Mesh);
        Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
        FixedPos = zeros(size(Mesh.Coordinates,1),1);
        FixedPos(Loc) = 1;
        NewMesh = jiggle(Mesh,FixedPos);   
      case 3
        Loc = get_BdEdges(Mesh);
        Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
        FixedPos = zeros(size(Mesh.Coordinates,1),1);
        FixedPos(Loc) = 1;
        NewMesh = smooth(Mesh,FixedPos);
    end
      
    % Assemble Stiffness matrix, load vector and incorporate BC
  
    A_LO = assemMat_QFE(NewMesh,@STIMA_Heat_QFE,P3O2(),SIGMA_HANDLE);
    L_LO = assemLoad_QFE(NewMesh,P3O2(),F_HANDLE);
    
    A_HI = assemMat_QFE(NewMesh,@STIMA_Heat_QFE,P7O6(),SIGMA_HANDLE);
    L_HI = assemLoad_QFE(NewMesh,P7O6(),F_HANDLE);
    
    % Incorporate Dirichlet and Neumann boundary data
         
    [U_LO,FreeDofs] = assemDir_QFE(NewMesh,-1,GD_HANDLE);
    L_LO = L_LO - A_LO*U_LO;
    
    [U_HI,FreeDofs] = assemDir_QFE(NewMesh,-1,GD_HANDLE);
    L_HI = L_HI - A_HI*U_HI;
    
    % Solve the linear system
  
    U_LO(FreeDofs) = A_LO(FreeDofs,FreeDofs)\L_LO(FreeDofs);
    U_HI(FreeDofs) = A_HI(FreeDofs,FreeDofs)\L_HI(FreeDofs);
       
    % Compute discretization error
    
    H1SErr_LO(i) = H1SErr_QFE(NewMesh,U_LO,P7O6(),U_EX);    
    H1SErr_HI(i) = H1SErr_QFE(NewMesh,U_HI,P7O6(),U_EX);
    h(i) = get_MeshWidth(Mesh);
    N(i) = size(NewMesh.Coordinates,1) + size(NewMesh.Edges,1);
    
  end
  
  % Plot out H1 semi-norm discretization error against h mesh width and add
  % slope triangles
  
  fig = figure('Name','Discretization errors');
  plot(h,H1SErr_LO,'r-', ...
       h,H1SErr_HI,'b-', ...
       h,H1SErr_LO,'k+', ...
       h,H1SErr_HI,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to H^1 semi-norm for quaratic finit elements}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Low order quadrature','High order quadrature','Location','NorthEast');
  p = polyfit(log(h),log(H1SErr_LO),1);
  add_Slope(gca,'North',p(1));
  p = polyfit(log(h),log(H1SErr_HI),1);
  add_Slope(gca,'East',p(1));
    
  % Plot out H1 semi-norm discretization error against number of dofs and
  % add slope triangles
  
  fig = figure('Name','Discretization errors');
  plot(N,H1SErr_LO,'r-', ...
       N,H1SErr_HI,'b-', ...
       N,H1SErr_LO,'k+', ...
       N,H1SErr_HI,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization errors with respect to H^1 semi-norm for quadratic finite elements}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Low order quadrature','High order quadrature','Location','NorthEast');
  p = polyfit(log(N),log(H1SErr_LO),1);
  add_Slope(gca,'North',p(1));
  p = polyfit(log(N),log(H1SErr_HI),1);
  add_Slope(gca,'SouthWest',p(1));
  
  % Clear memory
  
  clear all;
  