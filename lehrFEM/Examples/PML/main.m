% Convergence rates for piecewise quadratic and linear finite elements for
% the Poisson equation with Dirichlet boundary conditions on the square.
 
%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 7;                                                             % Number of red refinement steps
  
  % DATA 1
  F_HANDLE = @(x,varargin)-pi^2*...
      (2*cos(pi*x(:,1)).*cos(pi*x(:,2)).*(sin(pi*x(:,1))+sin(pi*x(:,2)))-...
       sin(pi*x(:,1)).*sin(pi*x(:,2)).*(x(:,1)+x(:,2)))-...
       sin(pi*x(:,1)).*sin(pi*x(:,2));               % Right hand side source term
  GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));             % Dirichlet boundary data
  U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));                % Exact solution for L2 norm
  EPS_Handle=@(x,varargin)[x(:,2) sin(pi*x(:,2)) sin(pi*x(:,1)) x(:,1)];
  
 % DATA 2
%     F_HANDLE = @(x,varargin)x(:,2).^2.*pi^2.*sin(pi*x(:,1))-...
%                             2*pi.*cos(pi*x(:,1)).*sin(pi*x(:,1))-...
%                             pi^2.*x(:,2).*cos(pi*x(:,1)).*cos(pi*x(:,2))-...
%                             pi*cos(pi*x(:,1)).*sin(pi*x(:,2));          % Right hand side source term
%     GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*x(:,2);                    % Dirichlet boundary data
%     U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*x(:,2);                       % Exact solution for L2 norm
%     EPS_Handle=@(x,varargin)[x(:,2) sin(pi*x(:,2)) sin(pi*x(:,1)) x(:,1)];

%   %DATA 3
%     F_HANDLE = @(x,varargin)-2*pi*cos(pi*x(:,1)).*sin(pi*x(:,2))-2*pi^2*...
%         (x(:,2).*cos(pi*x(:,1)).*cos(pi*x(:,2))-...
%          x(:,1).*sin(pi*x(:,1)).*sin(pi*x(:,2)))% Right hand side source term
%     GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));           % Dirichlet boundary data
%     U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));              % Exact solution for L2 norm
%     EPS_Handle=@(x,varargin)[x(:,1) x(:,2) x(:,2) x(:,1)];  
%    
    %DATA 4
%     F_HANDLE = @(x,varargin)2*i*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));     % Right hand side source term
%     GD_HANDLE = @(x,varargin)i*sin(pi*x(:,1)).*sin(pi*x(:,2));           % Dirichlet boundary data
%     U_EX_1 = @(x,varargin)i*sin(pi*x(:,1)).*sin(pi*x(:,2));              % Exact solution for L2 norm
%     EPS_Handle=@(x,varargin)[ones(size(x,1),1) zeros(size(x,1),1) zeros(size(x,1),1) ones(size(x,1),1)];  
%     MU_Handle=@(x,varagin)ones(size(x,1),1);

    %DATA 5
  F_HANDLE = @(x,varargin)-x(:,1)-x(:,2)-x(:,1).*x(:,2)    % Right hand side source term
  GD_HANDLE = @(x,varargin)-i*x(:,1).*x(:,2);           % Dirichlet boundary data
  U_EX_1 = @(x,varargin)-i*x(:,1).*x(:,2);
  EPS_Handle=@(x,varargin)i*[x(:,2)+x(:,1) zeros(size(x,1),1) zeros(size(x,1),1) x(:,2)];
  MU_Handle=@(x,varagin)i*ones(size(x,1),1);  
   
%  F_HANDLE=@(x,varargin)zeros(size(x,1),1);
%  GD_HANDLE=@DirData;
%  %GD_HANDLE=@(x,varargin)besselh(0,sqrt(x(:,1).^2+x(:,2).^2));
%  U_EX_1=@(x,varargin)besselh(0,sqrt(x(:,1).^2+x(:,2).^2));
%  EPS_Handle=@Epsilon;
%  MU_Handle=@Mu;
  JIG = 1;                                                               % Jiggle parameter
 
  % Initialize mesh
  QuadRule=Duffy(TProd(gauleg(0,1,10)));
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
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
      
    Mesh = refine_REG(Mesh);    
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
    %A_LFE = assemMat_LFE(NewMesh,@STIMA_LaplTensor_LFE,EPS_Handle,QuadRule);
    L_LFE = assemLoad_LFE(NewMesh,QuadRule,F_HANDLE);
    M_LFE = assemMat_LFE(NewMesh,@MASS_LFE);
    MW_LFE = assemMat_LFE(NewMesh,@MASS_WEIGHT_LFE,QuadRule,MU_Handle);
    UL_LFE = assemLoad_LFE(NewMesh,QuadRule,U_EX_1);
    % Incorporate Dirichlet and Neumann boundary data
%     u=M_LFE\L_LFE;
%     abs(u'*A_LFE*u)
    
    S=A_LFE-MW_LFE;
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(NewMesh,-1,GD_HANDLE);
    L_LFE = L_LFE - S*U_LFE;
      
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = S(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    
    
    plot_LFE(U_LFE,Mesh);
    %u=M_LFE\UL_LFE;
    colorbar;
    %Compute discretization error
    
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