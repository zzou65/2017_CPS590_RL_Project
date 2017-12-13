% Convergence rates for piecewise linear finite elements for
% Helmholtz equation using PML and Hankelfunction as Dirichlet Data
 
%   Copyright 2006-2006 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
%    NREFS =6;                                                         % Number of red refinement steps
%    F_HANDLE = @(x,varargin)-pi^2*...
%       (2*cos(pi*x(:,1)).*cos(pi*x(:,2)).*(sin(pi*x(:,1))+sin(pi*x(:,2)))-...
%        sin(pi*x(:,1)).*sin(pi*x(:,2)).*(x(:,1)+x(:,2)))-...
%        sin(pi*x(:,1)).*sin(pi*x(:,2));               % Right hand side source term
%    GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));           % Dirichlet boundary data
%    U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));              % Exact solution for L2 norm
%   EPS_Handle=@(x,varargin)[x(:,2) sin(pi*x(:,2)) sin(pi*x(:,1)) x(:,1)];
  
%   F_HANDLE = @(x,varargin)-x(:,1)-x(:,2)-x(:,1).*x(:,2)    % Right hand side source term
%   GD_HANDLE = @(x,varargin)x(:,1).*x(:,2);           % Dirichlet boundary data
%   U_EX_1 = @(x,varargin)x(:,1).*x(:,2);
%   EPS_Handle=@(x,varargin)[x(:,2)+x(:,1) zeros(size(x,1),1) zeros(size(x,1),1) x(:,2)];

%   F_HANDLE = @(x,varargin)x(:,2).^2.*pi^2.*sin(pi*x(:,1))-...
%                           2*pi.*cos(pi*x(:,1)).*sin(pi*x(:,1))-...
%                           pi^2.*x(:,2).*cos(pi*x(:,1)).*cos(pi*x(:,2))-...
%                           pi*cos(pi*x(:,1)).*sin(pi*x(:,2));          % Right hand side source term
%   GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*x(:,2);                      % Dirichlet boundary data
%   U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*x(:,2);                         % Exact solution for L2 norm
%   EPS_Handle=@(x,varargin)[x(:,2) sin(pi*x(:,2)) sin(pi*x(:,1)) x(:,1)];


%   F_HANDLE = @(x,varargin)-2*pi*cos(pi*x(:,1)).*sin(pi*x(:,2))-2*pi^2*...
%       (x(:,2).*cos(pi*x(:,1)).*cos(pi*x(:,2))-...
%        x(:,1).*sin(pi*x(:,1)).*sin(pi*x(:,2)))% Right hand side source term
%   GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));           % Dirichlet boundary data
%   U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));              % Exact solution for L2 norm
%   EPS_Handle=@(x,varargin)[x(:,1) x(:,2) x(:,2) x(:,1)];  
 
 F_HANDLE=@(x,varargin)zeros(size(x,1),1);
 GD_HANDLE=@DirData;
 %GD_HANDLE=@(x,varargin)besselh(0,sqrt(x(:,1).^2+x(:,2).^2));
 U_EX_1=@(x,varargin)besselh(0,sqrt(x(:,1).^2+x(:,2).^2));
 EPS_Handle=@Epsilon;
 MU_Handle=@Mu;
 slope=5*10.^[2 3];
 
  
%  F_HANDLE=@(x,varargin)zeros(size(x,1),1);
%  %GD_HANDLE=@DirData;
%  GD_HANDLE=@(x,varargin)exp(i*sqrt(x(:,1).^2+x(:,2).^2));
%  U_EX_1=@(x,varargin)exp(i*sqrt(x(:,1).^2+x(:,2).^2));
%  EPS_Handle=@Epsilon;
%  MU_Handle=@Mu;
    
 JIG =2;                                                               % Jiggle parameter
 NREFS =4;  
  
  % Initialize mesh
  QuadRule=Duffy(TProd(gauleg(0,1,10)));
  
  DHANDLE = inline('dist_diff(dist_rect(x,[-1 -1],2,2),dist_rect(x,[-0.25 -0.25],0.5,0.5))','x');
  
  Mesh = load_Mesh('pmlCoord.dat','pmlElem.dat');
  
  % define domain for error calculation [-a,a]^2
  a=0.75;
  Mesh.ElemFlag = max(abs([Mesh.Coordinates(Mesh.Elements(:,1),:) ...
                           Mesh.Coordinates(Mesh.Elements(:,2),:) ...
                           Mesh.Coordinates(Mesh.Elements(:,3),:)]),[],2)<=a;
                       
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  
  % calculate Boundary flags: -1 inner, 
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = (max(abs(Mesh.Coordinates(Mesh.Edges(Loc,1),:)),[],2)==0.25)-2;
  
  % Compute discretization error on a series of meshes
  nslope=size(slope,2); 
  h = zeros(nslope,NREFS);
  N_LFE = zeros(nslope,NREFS);
  L2_Error_LFE = zeros(nslope,NREFS);
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
    
    for s=1:nslope
    Sigma_Handle=@(x,x0) -slope(s)*(x-ones(size(x,1),1)*x0).^2;
    % Assemble Stiffness matrix, load vector and incorporate BC
    
    A_LFE = assemMat_LFE(NewMesh,@STIMA_hzm2,EPS_Handle,QuadRule,Sigma_Handle);
%   A_LFE = assemMat_LFE(NewMesh,@STIMA_LaplTensor_LFE,EPS_Handle,QuadRule);
    L_LFE = assemLoad_LFE(NewMesh,QuadRule,F_HANDLE);
    MW_LFE = assemMat_LFE(NewMesh,@MASS_WEIGHT_LFE,QuadRule,MU_Handle,Sigma_Handle);
%     M_LFE = assemMat_LFE(NewMesh,@MASS_LFE);
%     UL_LFE = assemLoad_LFE(NewMesh,QuadRule,U_EX_1);
    % Incorporate Dirichlet and Neumann boundary data
    S=A_LFE-MW_LFE;
    
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(NewMesh,[-1 -2],GD_HANDLE);
    L_LFE = L_LFE - S*U_LFE;
      
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = S(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    
%     plot_LFE(U_LFE,Mesh);
%     colorbar;
%     plot_LFE(M_LFE\UL_LFE,Mesh);
%     colorbar;
    % Compute discretization error
    
    L2_Error_LFE(s,i) = L2Err_LFE_splitDom(NewMesh,U_LFE,QuadRule,U_EX_1,1)
    %L2_Error_LFE(s,i) = L2Err_LFE(NewMesh,U_LFE,QuadRule,U_EX_1)
    N_LFE(s,i) = size(NewMesh.Coordinates,1);
 
    h(s,i) = get_MeshWidth(Mesh);
    end
  end
  
  % Plot out L2 discretization error against h mesh width and add slope
  % triangles
  
  fig = figure('Name','Discretization errors');
  for s=1:nslope
    plot(h(s,:),L2_Error_LFE(s,:));
    hold on;
  end
  hold off;
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Mesh width}');
  ylabel('{\bf Discretization error}');
 
  legend('Linear FE','Location','NorthEast');
  %p = polyfit(log(h),log(L2_Error_LFE),1);
  add_Slope(gca,'southwest',2);
    
    
  % Plot out L2 discretization error against number of dofs and add slope
  % triangles
  
  fig = figure('Name','Discretization errors');
 % set(gca,'ColorOrder',[1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1]);
  for s=1:nslope
      plot(N_LFE(s,:),L2_Error_LFE(s,:));
      hold on;
  end
  hold off;
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Dofs}');
  ylabel('{\bf Discretization error}');
  
  legend('Linear FE','Location','NorthEast');
 % p = polyfit(log(N_LFE(NREFS-3:NREFS)),log(L2_Error_LFE(NREFS-3:NREFS)),1);
  add_Slope(gca,'Southwest',-1);
   
  % Clear memory
  
  clear all;