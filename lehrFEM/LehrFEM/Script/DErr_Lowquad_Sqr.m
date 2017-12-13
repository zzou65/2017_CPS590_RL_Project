% Convergence rates for piecewise quadratic and linear finite elements for
% the Poisson equation with Dirichlet boundary conditions on the square.
% This script generates the following figures:
%   Q_lowquad_meshwidth.eps,  Q_lowquad_dofs.eps.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  INIT_NREFS = 1;                                                 % Number of initial mesh refinements
  NREFS = 6;                                                      % Number of red refinement steps
  F_HANDLE = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));        % Dirichlet boundary data
  U_EX = @(x,varargin)pi*[sin(pi*x(:,2)).*cos(pi*x(:,1)) ...      % Exact solution for H1 semi-norm
                          sin(pi*x(:,1)).*cos(pi*x(:,2))];
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Do initial mesh refinements
  
  for i = 1:INIT_NREFS
    Mesh = refine_REG(Mesh);
  end
  
  % Compute discretization error on a series of meshes
    
  h = zeros(1,NREFS);
  N = zeros(1,NREFS);
  H1S_Error = zeros(1,NREFS);
  H1S_Error_LO = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh);
    
    % Mesh preprocessing
    
    Loc = get_BdEdges(Mesh);
    Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
    FixedPos = zeros(size(Mesh.Coordinates,1),1);
    FixedPos(Loc) = 1;
    Mesh = jiggle(Mesh,FixedPos);   

    % Assemble stiffness matrix and load vector
  
    A = assemMat_QFE(Mesh,@STIMA_Lapl_QFE);
    L = assemLoad_QFE(Mesh,P7O6(),F_HANDLE);
    L_LO = assemLoad_QFE(Mesh,P3O2(),F_HANDLE);
    
    % Incorporate Dirichlet boundary conditions
         
    [U,FreeDofs] = assemDir_QFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    
    [U_LO,FreeDofs] = assemDir_QFE(Mesh,-1,GD_HANDLE);
    L_LO = L_LO - A*U_LO;
    
    % Solve the linear system
  
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    U_LO(FreeDofs) = A(FreeDofs,FreeDofs)\L_LO(FreeDofs);
    
    % Compute discretization error
    
    H1S_Error(i) = H1SErr_QFE(Mesh,U,P7O6(),U_EX);
    H1S_Error_LO(i) = H1SErr_QFE(Mesh,U_LO,P7O6(),U_EX);
    N(i) = size(Mesh.Coordinates,1) + size(Mesh.Edges,1);
    h(i) = get_MeshWidth(Mesh);
    
  end
  
  % Plot out H1 semi-norm discretization error against h mesh width and add
  % slope triangles
  
  fig = figure;
  plot(h,H1S_Error_LO,'r-', ...
       h,H1S_Error,'b-', ...
       h,H1S_Error_LO,'k+', ...
       h,H1S_Error,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to H^1 semi-norm for quaratic finit elements}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Low order quadrature','High order quadrature','Location','SouthWest');
  p = polyfit(log(h),log(H1S_Error_LO),1);
  add_Slope(gca,'NorthWest',p(1));
  p = polyfit(log(h),log(H1S_Error),1);
  add_Slope(gca,'South',p(1));
  
  print('-depsc','Q_lowquad_meshwidth.eps');
  close(fig);
  !gv Q_lowquad_meshwidth.eps &
    
  % Plot out H1 semi-norm discretization error against number of dofs and
  % add slope triangles
  
  fig = figure;
  plot(N,H1S_Error_LO,'r-', ...
       N,H1S_Error,'b-', ...
       N,H1S_Error_LO,'k+', ...
       N,H1S_Error,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization errors with respect to H^1 semi-norm for quadratic finite elements}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Low order quadrature','High order quadrature','Location','NorthEast');
  p = polyfit(log(N),log(H1S_Error_LO),1);
  add_Slope(gca,'North',p(1));
  p = polyfit(log(N),log(H1S_Error),1);
  add_Slope(gca,'SouthWest',p(1));
  
  print('-depsc','Q_lowquad_dofs.eps');
  close(fig);
  !gv Q_lowquad_dofs.eps &
  
  % Clear memory
  
  clear all;