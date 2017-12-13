% Convergence rates for piecewise quadratic and linear finite elements for
% the Poisson equation with Dirichlet boundary conditions on the square.
 
%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS =7;                                                         % Number of red refinement steps
%   F_HANDLE = @(x,varargin)-pi^2*...
%       (2*cos(pi*x(:,1)).*cos(pi*x(:,2)).*(sin(pi*x(:,1))+sin(pi*x(:,2)))-...
%        sin(pi*x(:,1)).*sin(pi*x(:,2)).*(x(:,1)+x(:,2)));               % Right hand side source term
%   GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));           % Dirichlet boundary data
%   U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));              % Exact solution for L2 norm
%   EPS_Handle=@(x,varargin)[x(:,2) sin(pi*x(:,2)) sin(pi*x(:,1)) x(:,1)];

  F_HANDLE = @(x,varargin)x(:,2).^2.*pi^2.*sin(pi*x(:,1))-...
                          2*pi.*cos(pi*x(:,1)).*sin(pi*x(:,1))-...
                          pi^2.*x(:,2).*cos(pi*x(:,1)).*cos(pi*x(:,2))-...
                          pi*cos(pi*x(:,1)).*sin(pi*x(:,2));          % Right hand side source term
  GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*x(:,2);                      % Dirichlet boundary data
  U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*x(:,2);                         % Exact solution for L2 norm
  EPS_Handle=@(x,varargin)[x(:,2) sin(pi*x(:,2)) sin(pi*x(:,1)) x(:,1)];


%   F_HANDLE = @(x,varargin)-2*pi*cos(pi*x(:,1)).*sin(pi*x(:,2))-2*pi^2*...
%       (x(:,2).*cos(pi*x(:,1)).*cos(pi*x(:,2))-...
%        x(:,1).*sin(pi*x(:,1)).*sin(pi*x(:,2)))% Right hand side source term
%   GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));           % Dirichlet boundary data
%   U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));              % Exact solution for L2 norm
%   EPS_Handle=@(x,varargin)[x(:,1) x(:,2) x(:,2) x(:,1)];  

  JIG = 2;                                                           % Jiggle parameter
 
  % Initialize mesh
  QuadRule=Duffy(TProd(gauleg(0,1,10)));
  
  DHANDLE = inline('dist_diff(dist_rect(x,[-2 -2],4,4),dist_circ(x,[0 0],0.5))','x')
  
  Mesh = load_Mesh('PMLSquareCoordc.dat','PMLSquareElemc.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Compute discretization error on a series of meshes
   
  h = zeros(1,NREFS);
  N_LFE = zeros(1,NREFS);
  L2_Error_LFE = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh,DHANDLE);    
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
    
    A_LFE = assemMat_LFE(NewMesh,@STIMA_hzm2,EPS_Handle,QuadRule);
    L_LFE = assemLoad_LFE(NewMesh,QuadRule,F_HANDLE);
    M_LFE = assemMat_LFE(NewMesh,@MASS_LFE);
    UL_LFE = assemLoad_LFE(NewMesh,QuadRule,U_EX_1);
    % Incorporate Dirichlet and Neumann boundary data
         
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(NewMesh,-1,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
      
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    
     plot_LFE(U_LFE,Mesh);
   %  plot_LFE(M_LFE\UL_LFE,Mesh);
     colorbar;
    % Compute discretization error
    
    L2_Error_LFE(i) = L2Err_LFE(NewMesh,U_LFE,QuadRule,U_EX_1)
    N_LFE(i) = size(NewMesh.Coordinates,1);
 
    h(i) = get_MeshWidth(Mesh);
    
  end
  
  % Plot out L2 discretization error against h mesh width and add slope
  % triangles
  
  fig = figure('Name','Discretization errors');
  plot(h,L2_Error_LFE,'b-');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Mesh width}');
  ylabel('{\bf Discretization error}');
  
  legend('Linear FE','Location','NorthEast');
  p = polyfit(log(h),log(L2_Error_LFE),1);
  add_Slope(gca,'North',p(1));
    
    
  % Plot out L2 discretization error against number of dofs and add slope
  % triangles
  
  fig = figure('Name','Discretization errors');
  plot(N_LFE,L2_Error_LFE,'b-');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Dofs}');
  ylabel('{\bf Discretization error}');
  
  legend('Linear FE','Location','NorthEast');
  p = polyfit(log(N_LFE(NREFS-3:NREFS)),log(L2_Error_LFE(NREFS-3:NREFS)),1);
  add_Slope(gca,'East',p(1));
   
  % Clear memory
  
  clear all;