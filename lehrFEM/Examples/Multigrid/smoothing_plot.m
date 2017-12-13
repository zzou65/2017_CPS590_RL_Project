function [] = smoothing_plot()
% plot errors for gauss-seidel smoother
%
%   This code plots the error distribution for a sample problem after 0, 1,
%   2 and 3 steps of Gauss-Seidel smoothing of a random initial guess.
%   Oscillation in the error are removed quickly, but there is little
%   apparent convergence to zero.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize constants
 
  ref = 4;                                       % refinements of initial mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  its = 3;                                       % number of smoothing iterations
  smoother = @gs_smooth;                         % multigrid smoother

  % generate mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i=1:ref
    Mesh = refine_REG(Mesh);
  end
  
  % assemble stiffness matrix and load vector
  
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  L = assemLoad_LFE(Mesh,P7O6(),f_handle);
  
  % incorporate boundary conditions
  
  [U_bd,FDofs] = assemDir_LFE(Mesh,-1,gd_handle);
  L = L - A*U_bd;
  A = A(FDofs,FDofs);
  L = L(FDofs);
  
  % calculate exact solution
  
  U_ex = U_bd;
  U_ex(FDofs) = A\L;
  
  % define oscillatory initial guess
  
  u0 = @(x) cos(0.5*pi*x(:,1)).*cos(0.5*pi*x(:,2))+sin(42*pi*x(:,1)).*sin(57*pi*x(:,2))+x(:,1).^2+x(:,2).^2;
  U = repmat(U_bd,1,its+1);
  U(FDofs,1) = u0(Mesh.Coordinates(FDofs,:));
  
  % initialize errors
  
  Err = zeros(size(U));
  Err(:,1) = U(:,1)-U_ex;
  
  % apply smoother
  
  for i=1:its
    U(FDofs,i+1) = smoother(U(FDofs,i),A,L);
    Err(:,i+1) = U(:,i+1)-U_ex;
  end
  
  % plot errors
  
  for i=1:its+1
    plot_LFE(Err(:,i),Mesh);
    set(gca,'CameraPosition',[1,2,1.5]);
    grid on;
    axis([-1,1,-1,1,-0.5,1])
    xlabel('\bf x');
    ylabel('\bf y');
    zlabel('\bf error');
    title(sprintf('%s Error After %.0f Smoothing Step(s)','\bf',i-1));
  end
  

return