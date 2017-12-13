%DERR_REG_REF-JIGGLE convergence rate of piecewise linear finite elements
% for the Poisson equation with Dirichlet boundary conditions on the square
% for meshes generated using refine_REG_jiggle with some jiggling.
%
% This code generates the following files:
%   DErr_reg_ref_jiggle.eps

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  INIT_NREFS = 1;                                                 % Number of initial mesh refinements
  NREFS = 6;                                                      % Number of red refinement steps
  JIGGLE = 0.25;                                                  % Amount of jiggling during refiniements
  F_HANDLE = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD_HANDLE = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));        % Dirichlet boundary data
  U_EX_1 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));           % Exact solution for L2 norm
  U_EX_2 = @(x,varargin)pi*[sin(pi*x(:,2)).*cos(pi*x(:,1)) ...    % Exact solution for H1 semi-norm
                            sin(pi*x(:,1)).*cos(pi*x(:,2))];
                          
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Do initial mesh refinements
  
  for i = 1:INIT_NREFS
    Mesh = refine_REG_jiggle(Mesh,JIGGLE);
  end
  
  % Compute discretization error on a series of meshes

  N = zeros(1,NREFS);
  L2_Error = zeros(1,NREFS);
  H1S_Error = zeros(1,NREFS);

  for i = 1:NREFS

    % Do red mesh refinement  

    Mesh = refine_REG_jiggle(Mesh,JIGGLE);
    
    % Assemble stiffness matrix and load vector

    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);

    % Incorporate Dirichlet boundary conditions

    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;

    % Solve the linear system

    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);

    % Compute discretization error

    L2_Error(i) = L2Err_LFE(Mesh,U,P7O6(),U_EX_1);
    H1S_Error(i) = H1SErr_LFE(Mesh,U,P7O6(),U_EX_2);
    N(i) = size(Mesh.Coordinates,1);
    
  end
  
  % Plot discretization errors against h and add slop triangles
  h = 1./sqrt(N);
  fig = figure;
  plot(h,L2_Error, h,H1S_Error);
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors for jiggled meshes}');
  xlabel('{\bf mesh width}');
  ylabel('{\bf discretization error}');
  legend('L_2 norm','H^1 seminorm','Location','NorthEast');
  p = polyfit(log(h),log(L2_Error),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(h),log(H1S_Error),1);
  add_Slope(gca,'NorthEast',p(1));
  
  print('-depsc','DErr_reg_ref_jiggle.eps');
  close(fig);
  !gv DErr_reg_ref_jiggle.eps &

% Clear memory
  
  clear all;
