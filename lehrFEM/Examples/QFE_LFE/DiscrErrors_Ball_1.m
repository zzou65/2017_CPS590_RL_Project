% Convergence rates for piecewise quadratic and linear finite elements for
% the Laplace equation with Dirichlet boundary conditions on the unit
% ball.
 
%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  DHANDLE = @dist_circ;                                    % Signed distance function
  C = [0 0];                                               % Center of the circle
  R = 1;                                                   % Radius of the circle
  NREFS = 7;                                               % Number of red refinement steps
  F_HANDLE = @(x,varargin)4*ones(size(x,1),1);             % Right hand side source term
  GD_HANDLE = @(x,varargin)zeros(size(x,1),1);             % Dirichlet boundary data 
  U_EX_1 = @(x,varargin)1-x(:,1).^2-x(:,2).^2;             % Exact solution for L2 norm
  U_EX_2 = @(x,varargin)-2*x;                              % Exact solution for H1 semi norm
  U_EX_3 = @(x,varargin)deal(1-x(:,1).^2-x(:,2).^2,-2*x);  % Exact solution for H1 norm
  JIG = 1;                                                 % Jiggle parameter
   
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Compute discretization error on a series of meshes
    
  h = zeros(1,NREFS);
  N_LFE = zeros(1,NREFS);
  N_QFE = zeros(1,NREFS);
  H1S_Error_LFE = zeros(1,NREFS);
  H1S_Error_QFE = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh,DHANDLE,C,R);    
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
  
    A_QFE = assemMat_QFE(NewMesh,@STIMA_Lapl_QFE);
    L_QFE = assemLoad_QFE(NewMesh,P7O6(),F_HANDLE);
    
    A_LFE = assemMat_LFE(NewMesh,@STIMA_Lapl_LFE);
    L_LFE = assemLoad_LFE(NewMesh,P7O6(),F_HANDLE);
    
    % Incorporate Dirichlet boundary data

    [U_LFE,FreeDofs_LFE] = assemDir_LFE(NewMesh,-1,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
    
    [U_QFE,FreeDofs_QFE] = assemDir_QFE(NewMesh,-1,GD_HANDLE);
    L_QFE = L_QFE - A_QFE*U_QFE;
    
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    U_QFE(FreeDofs_QFE) = A_QFE(FreeDofs_QFE,FreeDofs_QFE)\L_QFE(FreeDofs_QFE);
       
    % Compute discretization error
    
    H1S_Error_LFE(i) = H1SErr_LFE(NewMesh,U_LFE,P7O6(),U_EX_2);
    N_LFE(i) = size(NewMesh.Coordinates,1);
    
    H1S_Error_QFE(i) = H1SErr_QFE(NewMesh,U_QFE,P7O6(),U_EX_2);
    N_QFE(i) = size(NewMesh.Coordinates,1) + size(NewMesh.Edges,1);
    
    h(i) = get_MeshWidth(Mesh);
    
  end
  
  % Plot out H1 semi-norm discretization error against h mesh width and add
  % slope triangles
  
  fig = figure('Name','Discretization error');
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
  add_Slope(gca,'North',p(1));
  
  % Plot out H1 semi-norm discretization error against number of dofs and
  % add slope triangles
  
  fig = figure('Name','Discretization error');
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
  add_Slope(gca,'North',p(1));
 
  % Clear memory
  
  clear all;
  