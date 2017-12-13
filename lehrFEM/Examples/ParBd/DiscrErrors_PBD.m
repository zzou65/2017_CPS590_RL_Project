% Discretization errors for quadratic finite elements with parabolic
% boundary approximation.
 
%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 7;             % Number of inital mesh refinements
  DHANDLE = @dist_circ;  % Signed distance function
  C = [0 0];             % Center of the circle
  R = 1;                 % Radius of the circle
  F_HANDLE = @f;         % Right hand side source term
  GD_HANDLE = @g_D;      % Dirichlet boundary data
  UEX_1 = @uex_L2;       % Exact solution for L2 norm
  UEX_2 = @uex_H1S;      % Exact solution for H1 semi-norm
  UEX_3 = @uex_H1;       % Exact solution for H1 norm
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Initialize quadrature rule
  
  QuadRule = P7O6();
  
  h = zeros(1,NREFS);
  LInf_Error_LFE = zeros(1,NREFS);
  L2_Error_LFE = zeros(1,NREFS);
  H1S_Error_LFE = zeros(1,NREFS);
  H1_Error_LFE = zeros(1,NREFS);
  LInf_Error_PBD = zeros(1,NREFS);
  L2_Error_PBD = zeros(1,NREFS);
  H1S_Error_PBD = zeros(1,NREFS);
  H1_Error_PBD = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Refine the mesh  
      
    Mesh = refine_REG(Mesh,DHANDLE,C,R);
    
    % Add random noise to mesh 
    
    Loc = get_BdEdges(Mesh);
    Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
    FixedPos = zeros(size(Mesh.Coordinates,1),1);
    FixedPos(Loc) = 1;
    Mesh = jiggle(Mesh,FixedPos);
    
    % Assemble mass matrix and load vector
    
    A_QFE = assemMat_QFE(Mesh,@STIMA_Lapl_QFE);
    L_QFE = assemLoad_QFE(Mesh,QuadRule,F_HANDLE);
    
    Mesh = add_Edge2Elem(Mesh);
    Mesh = add_ParBd(Mesh,DHANDLE,C,R);
    A_PBD = assemMat_PBD(Mesh,@STIMA_Lapl_PBD,QuadRule);
    L_PBD = assemLoad_PBD(Mesh,QuadRule,F_HANDLE);
    
    % Incorporate Dirichlet boundary conditions
  
    [U_PBD,FreeDofs] = assemDir_PBD(Mesh,-1,GD_HANDLE);
    L_PBD = L_PBD - A_PBD*U_PBD;
  
    [U_QFE,FreeDofs] = assemDir_QFE(Mesh,-1,GD_HANDLE);
    L_QFE = L_QFE - A_QFE*U_QFE;
    
    % Solve the linear system
  
    U_PBD(FreeDofs) = A_PBD(FreeDofs,FreeDofs)\L_PBD(FreeDofs);
    U_QFE(FreeDofs) = A_QFE(FreeDofs,FreeDofs)\L_QFE(FreeDofs);
     
    % Compute discretization errors
    
    LInf_Error_QFE(i) = LInfErr_QFE(Mesh,U_QFE,UEX_1);
    L2_Error_QFE(i) = L2Err_QFE(Mesh,U_QFE,QuadRule,UEX_1);
    H1S_Error_QFE(i) = H1SErr_QFE(Mesh,U_QFE,QuadRule,UEX_2);
    H1_Error_QFE(i) = H1Err_QFE(Mesh,U_QFE,QuadRule,UEX_3);    
    
    LInf_Error_PBD(i) = LInfErr_PBD(Mesh,U_PBD,UEX_1);
    L2_Error_PBD(i) = L2Err_PBD(Mesh,U_PBD,QuadRule,UEX_1);
    H1S_Error_PBD(i) = H1SErr_PBD(Mesh,U_PBD,QuadRule,UEX_2);
    H1_Error_PBD(i) = H1Err_PBD(Mesh,U_PBD,QuadRule,UEX_3);
    
    h(i) = get_MeshWidth(Mesh);
    
  end
   
  % Discretization error for straight and curved elements with respect to
  % LInf norm
  
  fig = figure('Name','Discretization error');
  plot(h,LInf_Error_QFE,'r-', ...
       h,LInf_Error_PBD,'b-', ...
       h,LInf_Error_QFE,'k+', ...
       h,LInf_Error_PBD,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors measured in L^{\infty} norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('straight elements','curved elements','Location','NorthEast');
  p = polyfit(log(h),log(LInf_Error_QFE),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h),log(LInf_Error_PBD),1);
  add_Slope(gca,'North',p(1));
  
  % Discretization error for straight and curved elements with respect to
  % L2 norm
  
  fig = figure('Name','Discretization error');
  plot(h,L2_Error_QFE,'r-', ...
       h,L2_Error_PBD,'b-', ...
       h,L2_Error_QFE,'k+', ...
       h,L2_Error_PBD,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors measured in L^2 norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('straight elements','curved elements','Location','NorthEast');
  p = polyfit(log(h),log(L2_Error_QFE),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h),log(L2_Error_PBD),1);
  add_Slope(gca,'North',p(1));
  
  % Discretization error for straight and curved elements with respect to
  % H1 semi-norm
  
  fig = figure('Name','Discretization error');
  plot(h,H1S_Error_QFE,'r-', ...
       h,H1S_Error_PBD,'b-', ...
       h,H1S_Error_QFE,'k+', ...
       h,H1S_Error_PBD,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors measured in H^1 semi-norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('straight elements','curved elements','Location','NorthEast');
  p = polyfit(log(h),log(H1S_Error_QFE),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h),log(H1S_Error_PBD),1);
  add_Slope(gca,'North',p(1));
  
  % Discretization error for straight and curved elements with respect to
  % H1 norm
  
  fig = figure('Name','Discretization error');
  plot(h,H1_Error_QFE,'r-', ...
       h,H1_Error_PBD,'b-', ...
       h,H1_Error_QFE,'k+', ...
       h,H1_Error_PBD,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors measured in H^1 norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('straight elements','curved elements','Location','NorthEast');
  p = polyfit(log(h),log(H1_Error_QFE),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h),log(H1_Error_PBD),1);
  add_Slope(gca,'North',p(1));
  
  % Clear memory
  
  clear all;
  