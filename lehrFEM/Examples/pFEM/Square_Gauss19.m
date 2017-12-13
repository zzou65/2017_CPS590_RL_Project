% Shows the convergence rate (H1-seminorm and L2-norm) 
% of solving the Poisson equation
% -div(grad(u)) = f on a squared shape 
% with zero Dirichlet boundary condition
% where f(x,y) = 2 * pi^2 * sin(pi*x) * sin(pi*y)
% using pFEM and Gaussian quadrature of order 19.

% Copyright 2010 Roman Fuchs (2010 April 23rd)
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

clear;
  % Initialize constants
  
  PMAX = 20;                                                      % Maximum polynomial degree
  F = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));         % Right hand side source term
  GD = @(x,varargin)zeros(size(x,1),1);                           % Dirichlet boundary data
  GRAD_UEX = @(x,varargin)pi*[cos(pi*x(:,1)).*sin(pi*x(:,2)) ...  % Gradient of exact solution
                              sin(pi*x(:,1)).*cos(pi*x(:,2))]; 
  UEX = @(x,varargin) sin(pi*x(:,1)).*sin(pi*x(:,2));
  
  % Initialize mesh

  Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;  
  Mesh = add_Edge2Elem(Mesh);
  
  % Prepare mesh for longest edge bisection 
  
  nDofs = zeros(1,PMAX);
  H1S_error = zeros(1,PMAX);
  L2_error = zeros(1,PMAX);
  
  for p = 1:PMAX
      
    % Assign polynomial degrees and build dof maps
  
    EDofs = (p-1)*ones(size(Mesh.Edges,1),1);
    if(p > 2)
      CDofs = (p-1)*(p-2)/2*ones(size(Mesh.Elements,1),1);
    else
      CDofs = zeros(size(Mesh.Elements,1),1);
    end
    Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs);
  
    % Build shape functions and quadrature rules

    QuadRule_1D = gauleg(0,1,2*p);
    Shap_1D = shap_hp([QuadRule_1D.x zeros(size(QuadRule_1D.x))],p);
    QuadRule_2D = Duffy(TProd(QuadRule_1D));
    Shap_2D = shap_hp(QuadRule_2D.x,p);
  
    % Assemble global load vector and mass matrix
  
    A = assemMat_hp(Mesh,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D);
    L = assemLoad_hp(Mesh,Elem2Dof,QuadRule_2D,Shap_2D,F);
  
    % Incoporate Dirichlet boundary conditions
  
    [U,FreeDofs] = assemDir_hp(Mesh,Elem2Dof,-1,QuadRule_1D,Shap_1D,GD);
    L = L - A*U;
  
    % Solve the linear system
  
    if(p > 1)
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % Compute discretization errors
    
    nDofs(p) = size(U,1);
    
    QuadRule_2D = P73O19(); % use 73 Gaussian points of order 19
    Shap_2D = shap_hp(QuadRule_2D.x,p); % get the shape function values at the points
    H1S_error(p) = H1SErr_hp(Mesh,U,Elem2Dof,QuadRule_2D,Shap_2D,GRAD_UEX);
    L2_error(p) = L2Err_hp(Mesh,U,Elem2Dof,QuadRule_2D,Shap_2D,UEX);
  end
    
  % Generate figure
  
  fig1 = figure('Name','Convergence rates for p-FEM'); hold on;
  plot((1:PMAX).^2,H1S_error,'r-o');
  plot((1:PMAX).^2,L2_error,'b-sq');
  title('{\bf Discretization errors for p-FEM}')
  xlabel('{\bf p^2 }');
  ylabel('{\bf Discretization error}');
  legend('H^1 semi-norm', 'L^2 norm', 'Location', 'East');
  set(gca,'XScale','lin','YScale','log');
  hold off;
  
  
  fig2 = figure('Name','Convergence rates for p-FEM'); hold on;
  plot(nDofs, H1S_error, 'r-o');
  plot(nDofs, L2_error, 'b-sq');
  title('{\bf Discretization errors for p-FEM}')
  xlabel('{\bf N }');
  ylabel('{\bf Discretization error}');
  legend('H^1 semi-norm', 'L^2 norm', 'Location', 'East');
  set(gca,'XScale','lin','YScale','log');
  hold off;
  % Clear memory
  
  clear all;
  