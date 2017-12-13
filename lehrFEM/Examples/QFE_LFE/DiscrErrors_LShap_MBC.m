% Convergence rates for piecewise quadratic and linear finite elements for
% the Laplace equation with mixed boundary conditions on the L-shaped
% domain.
 
%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 7;                % Number of red refinement steps
  F_HANDLE = @f_LShap;      % Right hand side source term
  GD_HANDLE = @g_D_LShap;   % Dirichlet boundary data 
  GN_HANDLE = @g_N_LShap;   % Neumann boundary data   
  U_EX_1 = @uex_LShap_L2;   % Exact solution for L2 norm
  U_EX_2 = @uex_LShap_H1S;  % Exact solution for H1 semi norm
  U_EX_3 = @uex_LShap_H1;   % Exact solution for H1 norm
  JIG = 1;                  % Jiggle parameter
   
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = [-1 -2 -7 -7 -3 -4];
  
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
    Mesh = add_Edge2Elem(Mesh);
    
    % Mesh preprocessing
    
    switch (JIG)
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
    
    % Incorporate Dirichlet and Neumann boundary data
  
    L_LFE = assemNeu_LFE(NewMesh,-1:-1:-4,L_LFE,gauleg(0,1,4),GN_HANDLE);
    L_QFE = assemNeu_QFE(NewMesh,-1:-1:-4,L_QFE,gauleg(0,1,4),GN_HANDLE);
       
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(NewMesh,-7,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
    
    [U_QFE,FreeDofs_QFE] = assemDir_QFE(NewMesh,-7,GD_HANDLE);
    L_QFE = L_QFE - A_QFE*U_QFE;
    
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    U_QFE(FreeDofs_QFE) = A_QFE(FreeDofs_QFE,FreeDofs_QFE)\L_QFE(FreeDofs_QFE);
       
    % Compute discretization error
    
    LInf_Error_LFE(i) = LInfErr_LFE(NewMesh,U_LFE,U_EX_1);
    L2_Error_LFE(i) = L2Err_LFE(NewMesh,U_LFE,P7O6(),U_EX_1);
    H1S_Error_LFE(i) = H1SErr_LFE(NewMesh,U_LFE,P7O6(),U_EX_2);
    N_LFE(i) = size(NewMesh.Coordinates,1);
    
    LInf_Error_QFE(i) = LInfErr_QFE(NewMesh,U_QFE,U_EX_1);
    L2_Error_QFE(i) = L2Err_QFE(NewMesh,U_QFE,P7O6(),U_EX_1);
    H1S_Error_QFE(i) = H1SErr_QFE(NewMesh,U_QFE,P7O6(),U_EX_2);
    N_QFE(i) = size(NewMesh.Coordinates,1) + size(NewMesh.Edges,1);
    
    h(i) = get_MeshWidth(Mesh);
    
  end
  
  % Plot out LInf discretization error against h mesh width and add slope
  % triangles
  
  fig = figure('Name','Discretization errors');
  plot(h,LInf_Error_QFE,'r-', ...
       h,LInf_Error_LFE,'b-', ...
       h,LInf_Error_QFE,'k+', ...
       h,LInf_Error_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to L^{\infty} norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Quadratic FE','Linear FE','Location','NorthEast');
  p = polyfit(log(h),log(LInf_Error_QFE),1);
  add_Slope(gca,'North',p(1));
  p = polyfit(log(h),log(LInf_Error_LFE),1);
  add_Slope(gca,'SouthEast',p(1));
   
  % Plot out L2 discretization error against h mesh width and add slope
  % triangles
  
  fig = figure('Name','Discretization errors');
  plot(h,L2_Error_QFE,'r-', ...
       h,L2_Error_LFE,'b-', ...
       h,L2_Error_QFE,'k+', ...
       h,L2_Error_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
 
  legend('Quadratic FE','Linear FE','Location','NorthEast');
  p = polyfit(log(h),log(L2_Error_QFE),1);
  add_Slope(gca,'North',p(1));
  p = polyfit(log(h),log(L2_Error_LFE),1);
  add_Slope(gca,'East',p(1));
  
  % Plot out H1 semi-norm discretization error against h mesh width and add
  % slope triangles
  
  fig = figure('Name','Discretization errors');
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
  add_Slope(gca,'North',p(1));
  p = polyfit(log(h),log(H1S_Error_LFE),1);
  add_Slope(gca,'East',p(1));
  
  % Plot out L inf discretization error against number of dofs and add
  % slope triangles
  
  fig = figure('Name','Discretization errors');
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
  add_Slope(gca,'East',p(1));
  p = polyfit(log(N_LFE),log(LInf_Error_LFE),1);
  add_Slope(gca,'SouthWest',p(1));
  
  % Plot out L2 discretization error against number of dofs and add slope
  % triangles
  
  fig = figure('Name','Discretization errors');
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
  add_Slope(gca,'East',p(1));
  p = polyfit(log(N_LFE),log(L2_Error_LFE),1);
  add_Slope(gca,'SouthWest',p(1));
  
  % Plot out H1 semi-norm discretization error against number of dofs and
  % add slope triangles
  
  fig = figure('Name','Discretization errors');
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
  add_Slope(gca,'East',p(1));
  p = polyfit(log(N_LFE),log(H1S_Error_LFE),1);
  add_Slope(gca,'SouthWest',p(1));
  
  % Clear memory
  
  clear all;
  