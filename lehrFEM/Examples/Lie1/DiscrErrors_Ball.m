% Convergence upind Method for convection diffusion
 
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  epsilon=0.001;                                               % difusion parameter
  beta=1;                                                  % convection paramter
  DHANDLE = @dist_circ;                                    % Signed distance function
  C = [0 0];                                               % Center of the circle
  R = 1;                                                   % Radius of the circle
  NREFS =4;                                                % Number of red refinement steps
  MU_HANDLE=@(x,varargin)1;
  F_HANDLE = @(x,varargin) pi.*cos(pi.*x(:,1)).*sinh(pi.*x(:,2))+ ...
       pi.*sin(pi.*x(:,1)).*cosh(pi.*x(:,2));            % Right hand side source term
  GD_HANDLE =@(x,varargin)sin(pi*x(:,1)).*sinh(pi*x(:,2));  % Dirichlet boundary data 
  V_HANDLE=@(x,varargin)ones(size(x,1),2);
  U_EX = @(x,varargin) sin(pi * x(:,1)).*sinh(pi * x(:,2)); % Exact solution for L2 norm
  Grad_U_EX =@(x,varargin) [...
       pi.*cos(pi.*x(:,1)).*sinh(pi.*x(:,2)), ...
       pi.*sin(pi.*x(:,1)).*cosh(pi.*x(:,2))];             % Exact solution for H1 semi norm
   
  JIG=1;
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Compute discretization error on a series of meshes
    
  h = zeros(1,NREFS);
  N_LFE = zeros(1,NREFS);
  L2_Error = zeros(1,NREFS);
  H1S_Error = zeros(1,NREFS);
  %Mesh = refine_REG(Mesh,DHANDLE,C,R);    
  for i = 1:NREFS
      
    
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
    
    % Laplace 
    M = assemMat_W1F(Mesh,@MASS_W1F,MU_HANDLE, P7O6());
    TopGrad=assemMat_TopGrad(Mesh);
    TopRot=assemMat_TopRot(Mesh);
    ContrOne=assemMat_ContrOne(Mesh,V_HANDLE);
    MassZero=assemMat_MassZeroD(Mesh);
  
    A_fd =(epsilon*TopGrad'*M*TopGrad+beta*MassZero*ContrOne*TopGrad);
  
  % source term
  
  L = assemLoad_LFE(NewMesh,P7O6(),F_HANDLE);
  
  % Direchlet boundary
  
  [U_fd,FreeDofs] = assemDir_LFE(NewMesh,-1,GD_HANDLE);
  
  L_fd = L - A_fd*U_fd;
   
  % solving system
  
  U_fd(FreeDofs) = A_fd(FreeDofs,FreeDofs)\L_fd(FreeDofs);
  
  % Compute discretization error
  QuadRule = Duffy(TProd(gauleg(0,1,5)));
  L2_Error(i) = L2Err_LFE(NewMesh,U_fd,QuadRule,U_EX);
  H1S_Error(i) = H1SErr_LFE(NewMesh,U_fd, QuadRule, Grad_U_EX);
  N_LFE(i) = size(NewMesh.Elements,1);
  h(i) = get_MeshWidth(NewMesh);
    
    
  % Do red mesh refinement  
    
  Mesh = refine_REG(Mesh,DHANDLE,C,R);    
  Mesh = add_Edge2Elem(Mesh);
        
  end

  % Plot out H1 semi-norm discretization error against number of dofs and
  % add slope triangles
 
  H1_Error=sqrt(L2_Error.^2+epsilon*H1S_Error.^2);
  fig = figure('Name','Discretization error');
  plot(N_LFE,L2_Error,'ro--', ...
       N_LFE,H1S_Error,'b*:',N_LFE,H1_Error,'g*:'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf Dofs}');
  ylabel('{\bf Error}');
  
  legend('L^2-error u','H^1S-error u','H^1-error u','Location','NorthEast');
  p = polyfit(log(N_LFE(NREFS-3:NREFS)),log(L2_Error(NREFS-3:NREFS)),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(N_LFE(NREFS-3:NREFS)),log(H1S_Error(NREFS-3:NREFS)),1);
  add_Slope(gca,'East',p(1));
  p = polyfit(log(N_LFE(NREFS-3:NREFS)),log(H1_Error(NREFS-3:NREFS)),1);
  add_Slope(gca,'SouthWest',p(1));
  
  fig = figure('Name','Discretization error2');
  plot(h,L2_Error,'ro--', ...
       h,H1S_Error,'b*:',h,H1_Error,'g*:'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf Dofs}');
  ylabel('{\bf Error}');
  
  legend('L^2-error u','H^1S-error u','H^1-error u','Location','NorthEast');
  p = polyfit(log(h(NREFS-3:NREFS)),log(L2_Error(NREFS-3:NREFS)),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h(NREFS-3:NREFS)),log(H1S_Error(NREFS-3:NREFS)),1);
  add_Slope(gca,'East',p(1));
  p = polyfit(log(h(NREFS-3:NREFS)),log(H1_Error(NREFS-3:NREFS)),1);
  add_Slope(gca,'SouthWest',p(1));
  % Clear memory
  
  clear all;
  