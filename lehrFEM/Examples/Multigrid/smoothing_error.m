function [] = smoothing_error()
% convergence of smoothing iteration
%
%   This code plots the error of a Gauss-Seidel iteration in the L2 and H1
%   (semi-)norms for a sample FEM problem.  The L2 convergence is very
%   slow, but the H1 convergence is fast for the first few iterations and
%   then slows down.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize constants
 
  ref = 7;                                       % refinements of initial mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  smoother = @gs_smooth;                         % multigrid smoother
  its = 11;                                      % number of smoothing iterations
  
  err_l2 = nan(1,its);
  err_h1 = nan(1,its);
  
  
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
  
  % assemble matrices for scalar products
  
  H1 = A;
  L2 = assemMat_LFE(Mesh,@MASS_LFE);
  
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
  U = U_bd;
  U(FDofs) = u0(Mesh.Coordinates(FDofs,:));
  
  % calculate initial errors
  
  Err = U-U_ex;
  err_l2(1) = sqrt(Err'*L2*Err);
  err_h1(1) = sqrt(Err'*H1*Err);
  
  % iterate
  
  for i=2:its
    
    % do smoothing step
    
    U(FDofs) = smoother(U(FDofs),A,L);
    
    % calculate errors
    
    Err = U-U_ex;
    err_l2(i) = sqrt(Err'*L2*Err);
    err_h1(i) = sqrt(Err'*H1*Err);
    
  end
  
  % plot errors
  
  figure;
  semilogy(0:its-1,err_l2/err_l2(1),'-o',0:its-1,err_h1/err_h1(1),'-^');
  grid on;
  legend('L^2 Error','H_0^1 Error','Location','SouthWest');
  xlabel('\bf number of iterations');
  ylabel('\bf relative error');
  title('\bf Convergence of Gauss-Seidel Smoother');
  
return