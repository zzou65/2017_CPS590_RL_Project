% Convergence rates for piecewise quadratic and linear finite elements for
% the Laplace equation with Dirichlet boundary conditions on the unit
% ball.
 
%   Copyright 2005-2006 Patrick Meury & Kah Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  DHANDLE = @dist_circ;                                    % Signed distance function
  C = [0 0];                                               % Center of the circle
  R = 1;                                                   % Radius of the circle
  NREFS =6;                                               % Number of red refinement steps
  F_HANDLE = @f_Ball;                                      % Right hand side source term
  GD_HANDLE = @g_D_Ball;                                   % Dirichlet boundary data 
  U_EX_1 = @uex_Ball_L2;                                   % Exact solution for L2 norm
  U_EX_2 = @uex_Ball_H1S;                                  % Exact solution for H1 semi norm
  J_EX_1=@jex_Ball_H1div;
  J_EX_DIV_Handle = @(x,varargin)zeros(size(x,1),1);
  JIG =2;                                                 % Jiggle parameter
   
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
  E_TR=zeros(1,NREFS);
  L2_Error_LFE = zeros(1,NREFS);
  H1S_Error_LFE = zeros(1,NREFS);
  L2_Error_u = zeros(1,NREFS);
  H1S_Error_u = zeros(1,NREFS);
  L2_Error_j = zeros(1,NREFS);
  Hdiv_Error_j = zeros(1,NREFS);
  Mesh = refine_REG(Mesh,DHANDLE,C,R);    
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
    
    A = assemMat_FOSLS_TRP1(NewMesh,@STIMA_FOSLS_TRP1,P7O6());
    L = assemLoad_FOSLS_TRP1(NewMesh,P7O6(),F_HANDLE);
  
    A_LFE = assemMat_LFE(NewMesh,@STIMA_Lapl_LFE);
    L_LFE = assemLoad_LFE(NewMesh,P7O6(),F_HANDLE);
    
    % Incorporate Neumann boundary data
  
    % Incorporate Dirichlet boundary data
 
    [U,FreeDofs] = assemDir_FOSLS_TRP1(NewMesh,-1,GD_HANDLE);
    L = L - A*U;
  
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(NewMesh,-1,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
     
    % Solve the linear system
 
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    
       
    % Compute discretization error
    nEdges=size(NewMesh.Edges,1);
    nCoordinates = size(NewMesh.Coordinates,1);
    J=U(1:nEdges);
    U=U((nEdges+1):(nEdges+nCoordinates));
    QuadRule = Duffy(TProd(gauleg(0,1,5)));
    H1S_Error_LFE(i) = H1SErr_LFE(NewMesh,U_LFE,QuadRule,U_EX_2);  
    H1S_Error_u(i) = H1SErr_LFE(NewMesh,U,QuadRule,U_EX_2);
    L2_Error_u(i) = L2Err_LFE(NewMesh,U,QuadRule,U_EX_1);
    L2_Error_LFE(i) = L2Err_LFE(NewMesh,U_LFE,QuadRule,U_EX_1);
    L2_Error_j(i) = L2Err_W1F(NewMesh,J,QuadRule,J_EX_1);
    Hdiv_Error_j(i) = HCurlSErr_W1F(NewMesh,J,P7O6,J_EX_DIV_Handle);
    N_LFE(i) = size(NewMesh.Coordinates,1);
    E_TR(i) = size(NewMesh.Edges,1);
    h(i) = get_MeshWidth(NewMesh);
        
  end
 
  save Balldata H1S_Error_u H1S_Error_LFE L2_Error_u L2_Error_LFE L2_Error_j Hdiv_Error_j N_LFE E_TR h
  % Plot out H1 semi-norm discretization error against number of dofs and
  % add slope triangles
 
  fig = figure('Name','Discretization error');
  plot(N_LFE,L2_Error_u,'ro--', ...
       N_LFE,H1S_Error_u,'go--',...
       N_LFE,L2_Error_LFE,'b*:', ...
       N_LFE,H1S_Error_LFE,'c*:',...
       E_TR,L2_Error_j,'mo-',...
       E_TR,sqrt(Hdiv_Error_j.^2+L2_Error_j.^2),'y*-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf Dofs}');
  ylabel('{\bf Error}');
  
 legend('FOSLS: L^2-error u','FOSLS: H^1_{semi}-error u','LFE: L^2-error u','LFE: H^1_{semi}-error u','FOSLS: L^2-error j','FOSLS: Hdiv-error j','Location','NorthEast');
  p = polyfit(log(N_LFE),log(H1S_Error_u),1);
  add_Slope(gca,'NorthEast',p(1));
 
  % Clear memory
  
  clear all;
  