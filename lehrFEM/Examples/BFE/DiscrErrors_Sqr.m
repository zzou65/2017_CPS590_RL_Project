% Convergence rates for the Dirichlet problem to the Laplace equation using
% bilinear finite elements.
 
%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 7;                                                        % Number of red refinement steps
  F_HANDLE = @(x,varargin)2*pi^2*sin(pi*x(:,2)).*sin(pi*x(:,1));    % Right hand side source term
  GD_HANDLE = @(x,varargin)zeros(size(x,1),1);                      % Dirichlet boundary data  
  UEX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));              % Exact solution for L2 norm
  UEX_2 = @(x,varargin)deal(sin(pi*x(:,1)).*sin(pi*x(:,2)), ...     % Exact solution for H1 norm
                            pi*[cos(pi*x(:,1)).*sin(pi*x(:,2)) ...  
                                sin(pi*x(:,1)).*cos(pi*x(:,2))]);
              
  % Initialize quadrature rule
  
  QuadRule = TProd(gauleg(0,1,2));
                            
  % Initialize mesh
  
  Mesh.Coordinates = [0 0; 1 0; 1 1; 0 1];
  Mesh.Elements = [1 2 3 4];
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Compute discretization error on a series of meshes
    
  h = zeros(1,NREFS);
  L2Err = zeros(1,NREFS);
  H1Err = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh);    
            
    % Assemble Stiffness matrix, load vector and incorporate BC
  
    A = assemMat_BFE(Mesh,@STIMA_Lapl_BFE,QuadRule);
    L = assemLoad_BFE(Mesh,QuadRule,F_HANDLE);
       
    % Incorporate Dirichlet boundary data
    
    [U,FreeDofs] = assemDir_BFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;    
    
    % Solve the linear system
  
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);   
 
    % Compute discretization error
    
    L2Err(i) = L2Err_BFE(Mesh,U,QuadRule,UEX_1);
    H1Err(i) = H1Err_BFE(Mesh,U,QuadRule,UEX_2);
    
    h(i) = get_MeshWidth(Mesh);
    
  end
  
  % Plot out L2 and H1 discretization errors against h mesh width and add slope
  % triangles
  
  fig = figure;
  plot(h,L2Err,'r-', ...
       h,H1Err,'b-', ...
       h,L2Err,'k+', ...
       h,H1Err,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors for bilinear finite elements}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('L^2 norm','H^1 norm','Location','NorthEast');
  p = polyfit(log(h),log(L2Err),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h),log(H1Err),1);
  add_Slope(gca,'West',p(1));
      
  % Clear memory
  
  clear all;
  