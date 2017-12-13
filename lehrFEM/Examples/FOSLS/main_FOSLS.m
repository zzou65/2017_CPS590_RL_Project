% Run script for piecewise linear finite element solver.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
%   NREFS = 3;               % Number of red refinement steps
%   F_HANDLE = @f_LShap;     % Right hand side source term
%   GD_HANDLE = @g_D_LShap;  % Dirichlet boundary data
%   %GN_HANDLE = @g_N_LShap;  % Neumann boundary data
%  
%   % Initialize mesh
%   
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
%   
  NREFS =4;               % Number of red refinement steps
  C = [0 0];                    % Center of the circle
  R = 1;                        % Radius of the circle
  F_HANDLE=@(x,varargin)4*ones(size(x,1),1);        % righthandside
  MU_HANDLE=@(x,varargin)ones(size(x,1),1);        % righthandside
  GD_HANDLE=@(x,varargin)1-x(:,1).^2-x(:,2).^2;      % Direchletboundary Condition
  uex = @(x,varargin)1-x(:,1).^2-x(:,2).^2;          % Exact solution for L2 norm
  graduex = @(x,varargin)-2*x;                      % Exact solution for H1 semi norm
  jex_rot = @(x,varargin)[-2.*x(:,2),2.*x(:,1)];                      % Exact solution for H1 semi norm
  % Initialize mesh
   Mesh.Coordinates = [0 0; ...
                       1 0; ...
                       1  1; ... 
                      0  1];
  Mesh.Elements = [1 2 3; ...
                   1 3 4];
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  
   h = zeros(1,NREFS);
   N_LFE = zeros(1,NREFS);
   E_TR=zeros(1,NREFS);
   L2_Error_u = zeros(1,NREFS);
   H1S_Error_u = zeros(1,NREFS);
   L2_Error_j = zeros(1,NREFS);
   Hdiv_Error_j = zeros(1,NREFS);
    for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh);    
    Mesh = add_Edge2Elem(Mesh);
    
    % Mesh preprocessing
    
   
    NewMesh = Mesh;      
    % Assemble Stiffness matrix, load vector and incorporate BC
    
    A = assemMat_FOSLS_TRP1(Mesh,@STIMA_FOSLS_TRP1,P7O6());
    L = assemLoad_FOSLS_TRP1(Mesh,P7O6(),F_HANDLE);
    % Incorporate Neumann boundary data
 
    % Incorporate Dirichlet boundary data
 
    [U,FreeDofs] = assemDir_FOSLS_TRP1(Mesh,-1,GD_HANDLE);
    L = L - A*U;
  
    % Solve the linear system
 
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
    %Testing operator
%     Mj=assemMat_W1F(Mesh,@MASS_W1F, MU_HANDLE,P7O6());
%     Lj_rot=assemLoad_W1F(Mesh,P7O6(),jex_rot);
%     Lgradu=assemLoad_W1F(Mesh,P7O6(),graduex);
%     J_sol_rot=Mj\Lj_rot;
%     J_sol=Mj\Lgradu;
%     Mu=assemMat_LFE(Mesh,@MASS_LFE);
%     Lu=assemLoad_LFE(Mesh,P7O6(),uex);
%     U_sol=Mu\Lu;  
    
    % Compute discretization error
    nEdges=size(NewMesh.Edges,1);
    nCoordinates = size(NewMesh.Coordinates,1);
%     J_sol_rot'*A(1:nEdges,1:nEdges)* J_sol_rot
%     J_sol'*A(1:nEdges,1:nEdges)* J_sol-8/3

    %     U_sol'*A(nEdges+1:nEdges+nCoordinates,1:nEdges)*J_sol_rot-8/3
    %     J_sol'*A(1:nEdges,nEdges+1:nEdges+nCoordinates)*U_sol-8/3
    %    U_sol'*A(nEdges+1:nEdges+nCoordinates,nEdges+1:nEdges+nCoordinates)*U_sol-8/3
    J=U(1:nEdges);
    U=U((nEdges+1):(nEdges+nCoordinates));
      QuadRule = Duffy(TProd(gauleg(0,1,10)));
    H1S_Error_u(i) = H1SErr_LFE(NewMesh,U,QuadRule(),graduex);
    L2_Error_u(i) = L2Err_LFE(NewMesh,U,QuadRule(),uex);
    L2_Error_j(i) = L2Err_W1F(NewMesh,J,QuadRule(),jex_rot);
    Hdiv_Error_j(i) = HCurlSErr_W1F(Mesh,J,P7O6(),F_HANDLE);
      
    N_LFE(i) = size(NewMesh.Coordinates,1);
    E_TR(i) = size(NewMesh.Edges,1);
    h(i) = get_MeshWidth(NewMesh);
%     
%      plot_LFE(U,NewMesh);
%      colorbar;
%      plot_LFE(uex(NewMesh.Coordinates),NewMesh);
%      colorbar;
%      figure;
%      plot_W1F(J_sol_rot,NewMesh);
%      figure;
%      plot_W1F(J_sol,NewMesh);
     
   end
  
  % Plot out H1 semi-norm discretization error against h mesh width and add
  % slope triangles
  
  fig = figure('Name','Discretization error');
  plot(N_LFE,L2_Error_u,'b*-', ...
       N_LFE,H1S_Error_u,'k*-',...
       E_TR,L2_Error_j,'r*-',...
       E_TR,Hdiv_Error_j,'g*-');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to H^1 semi-norm}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
 legend('LFE: L^2-error','LFE: H^1-error','TR: L^2-error','TR: Hdiv-error','Location','NorthEast');
  p = polyfit(log(N_LFE),log(H1S_Error_u),1);
  add_Slope(gca,'North',p(1));
  
  
  clear all;
  