% Convergence rates for piecewise quadratic and linear finite elements for
% the Poisson equation with Dirichlet boundary conditions on the unit
% ball. This script generates the following .eps files:
%  B_meshwidth.eps,  B_dofs.eps.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  DHANDLE = @dist_circ;                            % Signed distance function
  C = [0 0];                                       % Center of the circle
  R = 1;                                           % Radius of the circle
  INIT_NREFS = 1;                                  % Number of initial mesh refinemenets
  NREFS = 6;                                       % Number of red refinement steps
  F_HANDLE = @(x,varargin)4*ones(size(x,1),1);     % Right hand side source term
  GD_HANDLE = @(x,varargin)1-x(:,1).^2-x(:,2).^2;  % Dirichlet boundary data 
  U_EX = @(x,varargin)-2*x;                        % Gradient of exact solution
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Do initial mesh refinements
  
  for i = 1:INIT_NREFS
    Mesh = refine_REG(Mesh,DHANDLE,C,R);
  end
  
  % Compute discretization error on a series of meshes
    
  h = zeros(1,NREFS);
  N_LFE = zeros(1,NREFS);
  N_QFE = zeros(1,NREFS);
  H1S_Error_LFE = zeros(1,NREFS);
  H1S_Error_QFE = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh,DHANDLE,C,R);    
    
    % Assemble Stiffness matrix, load vector and incorporate BC
  
    A_QFE = assemMat_QFE(Mesh,@STIMA_Lapl_QFE);
    L_QFE = assemLoad_QFE(Mesh,P7O6(),F_HANDLE);
    
    A_LFE = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L_LFE = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    
    % Incorporate Dirichlet boundary data

    [U_LFE,FreeDofs_LFE] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
    
    [U_QFE,FreeDofs_QFE] = assemDir_QFE(Mesh,-1,GD_HANDLE);
    L_QFE = L_QFE - A_QFE*U_QFE;
    
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    U_QFE(FreeDofs_QFE) = A_QFE(FreeDofs_QFE,FreeDofs_QFE)\L_QFE(FreeDofs_QFE);
     
    % Compute discretization error
    
    H1S_Error_LFE(i) = H1SErr_LFE(Mesh,U_LFE,P7O6(),U_EX);
    N_LFE(i) = size(Mesh.Coordinates,1);
    
    H1S_Error_QFE(i) = H1SErr_QFE(Mesh,U_QFE,P7O6(),U_EX);
    N_QFE(i) = size(Mesh.Coordinates,1) + size(Mesh.Edges,1);
    
    h(i) = get_MeshWidth(Mesh);
    
  end
  
  % Plot out H1 semi-norm discretization error against h mesh width and add
  % slope triangles
  
  fig = figure;
  plot(h,H1S_Error_QFE,'r-', ...
       h,H1S_Error_LFE,'b-', ...
       h,H1S_Error_QFE,'k+', ...
       h,H1S_Error_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to H^1 semi-norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Quadratic FE','Linear FE','Location','NorthEast');
  p = polyfit(log(h),log(H1S_Error_QFE),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h),log(H1S_Error_LFE),1);
  add_Slope(gca,'NorthEast',p(1));
  
  print('-depsc','B_meshwidth.eps');
  close(fig);
  !gv B_meshwidth.eps &
  
  % Plot out H1 semi-norm discretization error against number of dofs and
  % add slope triangles
  
  fig = figure;
  plot(N_QFE,H1S_Error_QFE,'r-', ...
       N_LFE,H1S_Error_LFE,'b-', ...
       N_QFE,H1S_Error_QFE,'k+', ...
       N_LFE,H1S_Error_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization errors with respect to H^1 semi-norm}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Quadratic FE','Linear FE','Location','NorthEast');
  p = polyfit(log(N_QFE),log(H1S_Error_QFE),1);
  add_Slope(gca,'SouthWest',p(1));
  p = polyfit(log(N_LFE),log(H1S_Error_LFE),1);
  add_Slope(gca,'NorthWest',p(1));
  
  print('-depsc','B_dofs.eps');
  close(fig);
  !gv B_dofs.eps &
  
  % Clear memory
  
  clear all;
  