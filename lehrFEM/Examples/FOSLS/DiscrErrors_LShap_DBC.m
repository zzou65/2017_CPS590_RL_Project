% Convergence rates for piecewise quadratic and linear finite elements for
% the Laplace equation with Dirichlet boundary conditions on the L-shaped
% domain.
 
%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS =6;                % Number of red refinement steps
  F_HANDLE = @f_LShap;      % Right hand side source term
  GD_HANDLE = @g_D_LShap;   % Dirichlet boundary data  
  U_EX_1 = @uex_LShap_L2;   % Exact solution for L2 norm
  U_EX_2 = @uex_LShap_H1S;  % Exact solution for H1 semi norm
  J_EX_1=@jex_LShap_H1div;
  ZERO_Handle = @(x,varargin)zeros(size(x,1),2);
  J_EX_DIV_Handle = @(x,varargin)zeros(size(x,1),1);
  JIG = 2;                  % Jiggle parameter
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Compute discretization error on a series of meshes
  
   E_TR=zeros(1,NREFS);
   h = zeros(1,NREFS);
   L2_Error_LFE = zeros(1,NREFS);
  H1S_Error_LFE = zeros(1,NREFS);
  N_LFE = zeros(1,NREFS);
  L2_Error_u = zeros(1,NREFS);
  H1S_Error_u = zeros(1,NREFS);
  L2_Error_j=zeros(1,NREFS);
  Hdiv_Error_j = zeros(1,NREFS);
  
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
    
    % Assemble stiffness matrix and load vector
 
  A = assemMat_FOSLS_TRP1(NewMesh,@STIMA_FOSLS_TRP1,P7O6());
  L = assemLoad_FOSLS_TRP1(NewMesh,P7O6(),F_HANDLE);
   
  A_LFE = assemMat_LFE(NewMesh,@STIMA_Lapl_LFE);
  L_LFE = assemLoad_LFE(NewMesh,P7O6(),F_HANDLE);
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
  
  
  QuadRule = Duffy(TProd(gauleg(0,1,10)));
  H1S_Error_LFE(i) = H1SErr_LFE(NewMesh,U_LFE,QuadRule,U_EX_2);  
  L2_Error_u(i) = L2Err_LFE(NewMesh,U,QuadRule(),U_EX_1);
  H1S_Error_u(i) = H1SErr_LFE(NewMesh,U,QuadRule(),U_EX_2);
  N_LFE(i) = size(NewMesh.Coordinates,1);
  L2_Error_j(i)= L2Err_W1F(NewMesh,J,QuadRule(),J_EX_1);      
  Hdiv_Error_j(i) = HCurlSErr_W1F(NewMesh,J,P7O6(),J_EX_DIV_Handle);
  L2_Error_LFE(i) = L2Err_LFE(NewMesh,U_LFE,QuadRule(),U_EX_1);
     
  E_TR(i) = size(NewMesh.Edges,1);
  h(i) = get_MeshWidth(NewMesh);
  end
  
  % Plot out H1 semi-norm discretization error against h mesh width and add
  % slope triangles
  
  
  % Plot out H1 semi-norm discretization error against number of dofs and
  % add slope triangles
   save LShapdata H1S_Error_u H1S_Error_LFE L2_Error_u L2_Error_LFE L2_Error_j Hdiv_Error_j N_LFE E_TR h
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
  
 legend('FOSLS: L^2-error u','FOSLS: H^1_{semi}-error u','LFE: L^2-error u','LFE: H_{semi}^1-error u','FOSLS: L^2-error j','FOSLS: Hdiv-error j','Location','NorthEast');
  p = polyfit(log(N_LFE),log(H1S_Error_u),1);
  add_Slope(gca,'NorthEast',p(1));
 
  % Clear memory
  
  clear all;
  
  