% Shows the L2-norm and H1-seminorm with respect to the polynomial degree
% p and the problem complexity N of the following problem:
%
% Solve the Poisson equation -div(grad(u)) = f on a L-Shape with
% u = g on the boundary (Dirichlet condition) using p-FEM.
% f is set to zero and g is the restriction of the exact solution u_ex
% to the boundary, where
% u_ex = r^(2/3) * sin(2/3 * phi),    in polar coordinates (r, phi).
%
% This problem is of interesting because grad(u) is discontinuous at
% the origin. 
% 
% This code is the same as /LehrFEM/hpFEM/main_4 with the following 
% extensions:
% - L2-norm is plotted as well
% - Gaussian quadrature of order 19 with 73 points is used to perform the
%   numerical integration
% 
% Copyright 2009 Roman Fuchs
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  clear;

  % Initialize constants  
  NREFS = 15;                                 % Number of mesh refinements  
  
  % Set up the differential equation and boundary data
  F        = @(x,varargin) zeros(size(x,1),1);% Right hand side source term
  UEX      = @uex_LShap;                      % exact solution on L-Shape
  GRAD_UEX = @grad_uex_LShap;                 % exact gradient of UEX
  GD       = @gD_LShap;                       % Dirichlet boundary data

  % Initialize mesh
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;  
  Mesh = add_Edge2Elem(Mesh);
  
  % Variables used for the plots at the end
  H1S_error = zeros(1,NREFS);
  L2_error  = zeros(1,NREFS);
  nDofs     = zeros(1,NREFS);
  
  for i = 1:NREFS 
    pmax=i; % do p-refinement
  
    EDofs = (pmax-1)*ones(size(Mesh.Edges,1),1);
    if(pmax > 2)
      CDofs = (pmax-1)*(pmax-2)/2*ones(size(Mesh.Elements,1),1);
    else
      CDofs = zeros(size(Mesh.Elements,1),1);
    end

    Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs);  
   
    % Build shape functions and quadrature rules
    QuadRule_1D = gauleg(0,1,2*pmax);
    Shap_1D = shap_hp([QuadRule_1D.x zeros(size(QuadRule_1D.x))],pmax);
    QuadRule_2D = Duffy(TProd(QuadRule_1D));
    Shap_2D = shap_hp(QuadRule_2D.x,pmax);
  
    % Assemble global load vector and mass matrix
    A = assemMat_hp(Mesh,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D);
    L = assemLoad_hp(Mesh,Elem2Dof,QuadRule_2D,Shap_2D,F);
  
    % Incoporate Dirichlet boundary conditions
    [U,FreeDofs] = assemDir_hp(Mesh,Elem2Dof,-1,QuadRule_1D,Shap_1D,GD);
    L = L - A*U;
  
    % Solve the linear system
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        
    % Compute discretization errors
    nDofs(i) = size(U,1);
    
    % Numerical integration with Gaussian quadrature of order 19
    QuadRule_2D = P73O19(); 
    Shap_2D = shap_hp(QuadRule_2D.x,i);
    H1S_error(i) = H1SErr_hp(Mesh,U,Elem2Dof,QuadRule_2D,Shap_2D,GRAD_UEX);
    L2_error(i)  = L2Err_hp( Mesh,U,Elem2Dof,QuadRule_2D,Shap_2D,UEX);
  end
  
  
  % Generate figures:
  % fig1: Norm discretization error w.r.t. polynomial degree (p)
  fig1 = figure('Name', 'Convergence rates for p-FEM'); hold on;
  plot(1:pmax, H1S_error, 'r-o');
  plot(1:pmax, L2_error, 'b-sq');
  title('{\bf Discretization errors for p-FEM}')
  xlabel('{\bf p}');
  ylabel('{\bf Norm of Discretization error}');
  legend('H^1 semi-norm', 'L^2 norm', 'Location', 'Northeast');
  hold off;
  set(gca, 'XScale', 'log', 'YScale', 'log');
  
  % fig2: Norm discretization error w.r.t. number of degrees of freedom (N)
  fig2 = figure('Name', 'Convergence rates for p-FEM'); hold on;
  plot(nDofs, H1S_error, 'r-o');
  plot(nDofs, L2_error, 'b-sq');
  title('{\bf Discretization errors for p-FEM}')
  xlabel('{\bf N }');
  ylabel('{\bf Norm of Discretization error}');
  legend('H^1 semi-norm', 'L^2 norm', 'Location', 'Northeast');
  hold off;
  set(gca,'XScale','log','YScale','log');
  
  % Clear memory
  clear all;
  
