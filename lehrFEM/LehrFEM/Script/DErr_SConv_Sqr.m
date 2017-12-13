% Convergence rates for piecewise quadratic finite elements for the Poisson
% equation with Dirichlet boundary conditions on structured grid of the
% unit square (superconvergence). This script generates the following
% figures:
%  Q_reg_meshwidth.eps,  Q_reg_dofs.eps.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  INIT_NREFS = 1;                                                 % Number of initial mesh refinements
  NREFS = 6;                                                      % Number of red refinement steps
  F_HANDLE = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));        % Dirichlet boundary data
  U_EX = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));             % Exact solution
  
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
  N_LFE = zeros(1,NREFS);
  N_QFE = zeros(1,NREFS);
  LInf_Error_LFE = zeros(1,NREFS);
  LInf_Error_QFE = zeros(1,NREFS);
  for i = 1:NREFS
  
    % Do red refinement
    
    Mesh = refine_REG(Mesh);
    
    % Assemble stiffness matrix and load vector
    
    A_LFE = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L_LFE = assemLoad_LFE(Mesh,P3O2(),F_HANDLE);
    
    A_QFE = assemMat_QFE(Mesh,@STIMA_Lapl_QFE);
    L_QFE = assemLoad_QFE(Mesh,P3O2(),F_HANDLE);
    
    % Incorporate Dirichlet boundary condition
    
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
    
    [U_QFE,FreeDofs_QFE] = assemDir_QFE(Mesh,-1,GD_HANDLE);
    L_QFE = L_QFE - A_QFE*L_QFE;
    
    % Solve the linear system
    
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    U_QFE(FreeDofs_QFE) = A_QFE(FreeDofs_QFE,FreeDofs_QFE)\L_QFE(FreeDofs_QFE);
    
    % Compute discretization error
    
    LInf_Error_LFE(i) = LInfErr_LFE(Mesh,U_LFE,U_EX);
    N_LFE(i) = size(Mesh.Coordinates,1);
    LInf_Error_QFE(i) = LInfErr_QFE(Mesh,U_QFE,U_EX);
    N_QFE(i) = size(Mesh.Coordinates,1)+size(Mesh.Edges,1);
    
    h(i) = get_MeshWidth(Mesh);
    
  end
  
  % Plot out H1 semi-norm discretization error with respect to h and add
  % slope triangles
  
  fig = figure;
  plot(h,LInf_Error_LFE,'r-', ...
       h,LInf_Error_QFE,'b-', ...
       h,LInf_Error_LFE,'k+', ...
       h,LInf_Error_QFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization error with respect to L^{\infty} norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Linear FE','Quadratic FE','Location','NorthWest');
  p = polyfit(log(h),log(LInf_Error_LFE),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h),log(LInf_Error_QFE),1);
  add_Slope(gca,'NorthWest',p(1));
  
  print('-depsc','Q_reg_meshwidth.eps');
  close(fig);
  !gv Q_reg_meshwidth.eps &
  
  % Plot out H1 semi-norm discretization error with respect to N and add
  % slope triangles
  
  fig = figure;
  plot(N_LFE,LInf_Error_LFE,'r-', ...
       N_QFE,LInf_Error_QFE,'b-', ...
       N_LFE,LInf_Error_LFE,'k+', ...
       N_QFE,LInf_Error_QFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization error with respect to L^{\infty} norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Linear FE','Quadratic FE','Location','NorthWest');
  p = polyfit(log(N_LFE),log(LInf_Error_LFE),1);
  add_Slope(gca,'NorthEast',p(1));
  p = polyfit(log(N_QFE),log(LInf_Error_QFE),1);
  add_Slope(gca,'SouthWest',p(1));
  
  print('-depsc','Q_reg_dofs.eps');
  close(fig);
  !gv Q_reg_dofs.eps &
  
  % Clear memory
  
  clear all;
  