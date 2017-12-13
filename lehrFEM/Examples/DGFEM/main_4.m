% Run script for discontinuous Galerkin finite element solver.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
  NREFS = 6;                                % Number of red refinement steps
  UEX = @(x,varargin)x(:,1).^2+x(:,2).^2;   % Exact solution
  GRAD_UEX = @(x,varargin)2*x;              % Gradient of exact solution
  G = @(x,varargin)-4*ones(size(x,1),1);    % Right hand side load data
  UD = @(x,varargin)x(:,1).^2+x(:,2).^2;    % Dirichlet boundary data
  S = 1;                                    % Symmetric (+1) or antisymmetric (-1) discretization
  SIGMA = @(P0,P1,varargin)10/norm(P1-P0);  % Edge weight function
  
  % Initialize mesh
  
  Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  
  h = zeros(1,NREFS);
  L2Err = zeros(1,NREFS);
  H1SErr = zeros(1,NREFS);
  QuadRule_1D = gauleg(0,1,2);
  QuadRule_2D = P3O3();
  for i = 1:NREFS
  
    Mesh = refine_REG(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Mesh = add_DGData(Mesh);
    
    % Assemble matrices and load vectors (discontinuous Lagrangian elements)
   
    [I1,J1,Avol] = assemMat_Vol_DG(Mesh,@STIMA_Lapl_Vol_DGLFE);
  
    [I2,J2,Jinn] = assemMat_Inn_DG(Mesh,@STIMA_InnPen_DGLFE,SIGMA);
    [I2,J2,Ainn] = assemMat_Inn_DG(Mesh,@STIMA_Inn_DGLFE,S);
    
    [I3,J3,Jbnd] = assemMat_Bnd_DG(Mesh,@STIMA_BndPen_DGLFE,SIGMA);  
    [I3,J3,Abnd] = assemMat_Bnd_DG(Mesh,@STIMA_Bnd_DGLFE,S);
  
    Lvol = assemLoad_Vol_DG(Mesh,@LOAD_Vol_DGLFE,QuadRule_2D,G);
    Lbnd = assemLoad_Bnd_DG(Mesh,@LOAD_Bnd_DGLFE,QuadRule_1D,S,SIGMA,UD);
    
    % Create system matrix
       
    A = sparse([J1; J2; J3; J2; J3], ...
               [I1; I2; I3; I2; I3], ...
               [Avol; Ainn; Abnd; Jinn; Jbnd]);
    L = Lvol + Lbnd;
         
    % Solve the linear system
  
    U = A\L;
    
    % Compute discretization error
    
    h(i) = get_MeshWidth(Mesh);
    L2Err(i) = L2Err_DGLFE(Mesh,U,QuadRule_2D,UEX);
    H1SErr(i) = H1SErr_DGLFE(Mesh,U,QuadRule_2D,GRAD_UEX);
    
  end
  
  % Generate figure
  
  fig = figure('Name','Discretization error');
  plot(h,L2Err,'r-', ...
       h,H1SErr,'b-', ...
       h,L2Err,'k+', ...
       h,H1SErr,'k+');
  title('{\bf Discretization error}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  set(gca,'XScale','log','XDir','reverse','YScale','log');
  
  legend('L^2 norm','H^1 semi-norm','Location','SouthWest');
  
  p = polyfit(log(h),log(L2Err),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(h),log(H1SErr),1);
  add_Slope(gca,'NorthWest',p(1));
  
  % Clear memory
  
  clear all;
  