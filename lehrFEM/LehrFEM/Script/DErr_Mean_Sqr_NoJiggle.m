% Convergence rates for piecewise quadratic and linear finite elements for
% the Poisson equation with Dirichlet boundary conditions on the square.
% This script generates the following figures:
%   Q_Mean_meshwidth_scvg.eps,     Q_Mean_meshwidth_scvg.eps,

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  INIT_NREFS = 1;                                                 % Number of initial mesh refinements
  NREFS = 6;                                                      % Number of red refinement steps
  F_HANDLE = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));        % Dirichlet boundary data          
  FEX = 4/pi^2;                                                   % Exact value of functional
  
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
  MeanErr_LFE = zeros(1,NREFS);
  MeanErr_QFE = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh);
    Mesh2 = Mesh;
    
    
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
    
    MeanErr_LFE(i) = abs(FEX - mean_LFE(U_LFE,Mesh2));
    N_LFE(i) = size(Mesh2.Coordinates,1);
    
    MeanErr_QFE(i) = abs(FEX - mean_QFE(U_QFE,Mesh2));
    N_QFE(i) = size(Mesh2.Coordinates,1) + size(Mesh2.Edges,1);
    
    h(i)  = get_MeshWidth(Mesh2);
    
  end  
  
  % Plot out mean error against h mesh width and add slope triangles
  
  fig = figure;
  plot(h,MeanErr_QFE,'r-', ...
       h,MeanErr_LFE,'b-', ...
       h,MeanErr_QFE,'k+', ...
       h,MeanErr_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Errors between the mean and exact solution}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Mean error [log]}');
  
  legend('Quadratic FE','Linear FE','Location','SouthWest');
  p = polyfit(log(h),log(MeanErr_QFE),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(h),log(MeanErr_LFE),1);
  add_Slope(gca,'NorthWest',p(1));
  
  print('-depsc','Q_Mean_meshwidth_scvg.eps');
  close(fig);
  !gv Q_Mean_meshwidth_scvg.eps &
  
  % Plot out mean error against N number of dofs and add slope triangles
  
  fig = figure;
  plot(N_QFE,MeanErr_QFE,'r-', ...
       N_LFE,MeanErr_LFE,'b-', ...
       N_QFE,MeanErr_QFE,'k+', ...
       N_LFE,MeanErr_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Error between mean and exact solution}');
  xlabel('{\bf Number of dofs [log]}');
  ylabel('{\bf Mean error [log]}');
  legend('Quadratic FE','Linear FE','Location','NorthEast');
  
  p = polyfit(log(N_QFE),log(MeanErr_QFE),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(N_LFE),log(MeanErr_LFE),1);
  add_Slope(gca,'North',p(1));
  
  print('-depsc','Q_Mean_dofs_scvg.eps');
  close(fig);
  !gv Q_Mean_dofs_scvg.eps & 
 
  % Clear memory
  
  clear all;
  