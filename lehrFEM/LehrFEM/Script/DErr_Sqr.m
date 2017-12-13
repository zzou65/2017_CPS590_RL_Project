% Convergence rates for piecewise quadratic and linear finite elements for
% the Poisson equation with Dirichlet boundary conditions on the square.
% This script generates the following figures:
%   Q_Linf_meshwidth.eps,     Q_Linf_dofs.eps,
%   Q_L2_meshwidth.eps,       Q_L2_dofs.eps,
%   Q_H1_meshwidth.eps,       Q_H1_dofs.eps.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  INIT_NREFS = 1;                                                 % Number of initial mesh refinements
  NREFS = 6;                                                      % Number of red refinement steps
  F_HANDLE = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));        % Dirichlet boundary data
  U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));           % Exact solution for L2 norm
  U_EX_2 = @(x,varargin)pi*[sin(pi*x(:,2)).*cos(pi*x(:,1)) ...    % Exact solution for H1 semi-norm
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
  N_LFE = zeros(1,NREFS);
  N_QFE = zeros(1,NREFS);
  LInf_Error_LFE = zeros(1,NREFS);
  LInf_Error_QFE = zeros(1,NREFS);
  L2_Error_LFE = zeros(1,NREFS);
  L2_Error_QFE = zeros(1,NREFS);
  H1S_Error_LFE = zeros(1,NREFS);
  H1S_Error_QFE = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh);
    Mesh2 = Mesh;
    
    % Mesh preprocessing
    
    Loc = get_BdEdges(Mesh2);
    Loc = unique([Mesh2.Edges(Loc,1); Mesh2.Edges(Loc,2)]);
    FixedPos = zeros(size(Mesh2.Coordinates,1),1);
    FixedPos(Loc) = 1;
    Mesh2 = jiggle(Mesh2,FixedPos);   

% Plot of perturbed mesh
     
     if (i <= 2)
       plot_Mesh(Mesh2,'as');
%       title(sprintf('Mesh on level %d',i),'fontsize',14);
       print('-depsc',sprintf('Q_mesh%d.eps',i));
     end

    % Assemble stiffness matrix and load vector
  
    A_QFE = assemMat_QFE(Mesh2,@STIMA_Lapl_QFE);
    L_QFE = assemLoad_QFE(Mesh2,P7O6(),F_HANDLE);
   
    A_LFE = assemMat_LFE(Mesh2,@STIMA_Lapl_LFE);
    L_LFE = assemLoad_LFE(Mesh2,P7O6(),F_HANDLE);
        
    % Incorporate Dirichlet boundary conditions
         
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(Mesh2,-1,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
    
    [U_QFE,FreeDofs_QFE] = assemDir_QFE(Mesh2,-1,GD_HANDLE);
    L_QFE = L_QFE - A_QFE*U_QFE;
    
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    U_QFE(FreeDofs_QFE) = A_QFE(FreeDofs_QFE,FreeDofs_QFE)\L_QFE(FreeDofs_QFE);
    
    % Compute discretization error
    
    LInf_Error_LFE(i) = LInfErr_LFE(Mesh2,U_LFE,U_EX_1);
    L2_Error_LFE(i) = L2Err_LFE(Mesh2,U_LFE,P7O6(),U_EX_1);
    H1S_Error_LFE(i) = H1SErr_LFE(Mesh2,U_LFE,P7O6(),U_EX_2);
    N_LFE(i) = size(Mesh2.Coordinates,1);
    
    LInf_Error_QFE(i) = LInfErr_QFE(Mesh2,U_QFE,U_EX_1);
    L2_Error_QFE(i) = L2Err_QFE(Mesh2,U_QFE,P7O6(),U_EX_1);
    H1S_Error_QFE(i) = H1SErr_QFE(Mesh2,U_QFE,P7O6(),U_EX_2);
    N_QFE(i) = size(Mesh2.Coordinates,1) + size(Mesh2.Edges,1);
    
    h(i) = get_MeshWidth(Mesh2);
    
  end

  % Plot out L inf discretization error against h mesh width and add slope
  % triangles
  
  fig = figure;
  plot(h,LInf_Error_QFE,'r-', ...
       h,LInf_Error_LFE,'b-', ...
       h,LInf_Error_QFE,'k+', ...
       h,LInf_Error_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to L^{\infty} norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Quadratic FE','Linear FE','Location','SouthWest');
  p = polyfit(log(h),log(LInf_Error_QFE),1);
  add_Slope(gca,'East',p(1));
  p = polyfit(log(h),log(LInf_Error_LFE),1);
  add_Slope(gca,'NorthWest',p(1));
  
  print('-depsc','Q_Linf_meshwidth.eps');
  close(fig);
  !gv Q_Linf_meshwidth.eps &
  
  % Plot out L2 discretization error against h mesh windth and add slope
  % triangles
  
  fig = figure;
  plot(h,L2_Error_QFE,'r-', ...
       h,L2_Error_LFE,'b-', ...
       h,L2_Error_QFE,'k+', ...
       h,L2_Error_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Quadratic FE','Linear FE','Location','SouthWest');
  p = polyfit(log(h),log(L2_Error_QFE),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(h),log(L2_Error_LFE),1);
  add_Slope(gca,'NorthWest',p(1));
  
  print('-depsc','Q_L2_meshwidth.eps');
  close(fig);
  !gv Q_L2_meshwidth.eps &
  
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
  
  legend('Quadratic FE','Linear FE','Location','SouthWest');
  p = polyfit(log(h),log(H1S_Error_QFE),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(h),log(H1S_Error_LFE),1);
  add_Slope(gca,'NorthWest',p(1));
  
  print('-depsc','Q_H1_meshwidth.eps');
  close(fig);
  !gv Q_H1_meshwidth.eps &
  
  % Plot out L inf discretization error against number of dofs and add
  % slope triangles
  
  fig = figure;
  plot(N_QFE,LInf_Error_QFE,'r-', ...
       N_LFE,LInf_Error_LFE,'b-', ...
       N_QFE,LInf_Error_QFE,'k+', ...
       N_LFE,LInf_Error_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization errors with respect to L^{\infty} norm}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  legend('Quadratic FE','Linear FE','Location','NorthEast');
    
  p = polyfit(log(N_QFE),log(LInf_Error_QFE),1);
  add_Slope(gca,'South',p(1))
  p = polyfit(log(N_LFE),log(LInf_Error_LFE),1);
  add_Slope(gca,'East',p(1));
  
  print('-depsc','Q_Linf_dofs.eps');
  close(fig);
  !gv Q_Linf_dofs.eps &
  
  % Plot out L2 discretization error against number of dofs and add slope
  % triangles
  
  fig = figure;
  plot(N_QFE,L2_Error_QFE,'r-', ...
       N_LFE,L2_Error_LFE,'b-', ...
       N_QFE,L2_Error_QFE,'k+', ...
       N_LFE,L2_Error_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  legend('Quadratic FE','Linear FE','Location','NorthEast');
  
  p = polyfit(log(N_QFE),log(L2_Error_QFE),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(N_LFE),log(L2_Error_LFE),1);
  add_Slope(gca,'East',p(1));
  
  print('-depsc','Q_L2_dofs.eps');
  close(fig);
  !gv Q_L2_dofs.eps &
  
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
  add_Slope(gca,'South',p(1));
  p = polyfit(log(N_LFE),log(H1S_Error_LFE),1);
  add_Slope(gca,'North',p(1));
  
  print('-depsc','Q_H1_dofs.eps');
  close(fig);
  !gv Q_H1_dofs.eps &
  
  % Clear memory
  
  clear all;
  