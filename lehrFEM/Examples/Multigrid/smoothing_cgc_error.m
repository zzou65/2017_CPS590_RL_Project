function [] = smoothing_cgc_error()
% convergence of smoothing with coarse grid correction
%
%   This code plots the error in the L2 norm and H1 seminorm after
%   presmoothing, a subsequent coarse grid correction and postsmoothing for
%   a sample poblem.  The smoothing decreases the H1 seminorm much more
%   than the L2 norm; fittingly, the coarse grid correction greatly
%   decreases the L2 norm and has a much smaller effect on the H1 seminorm.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize constants
 
  ref = 7;                                       % refinements of initial mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  smoother = @gs_smooth;                         % multigrid smoother
  m = 2;                                         % number of smoothing iterations
  
  err_l2 = nan(1,3*m+3);
  err_h1 = nan(1,3*m+3);
  
  % generate mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i=1:ref
    cMesh = Mesh;
    Mesh = refine_REG(cMesh);
  end
  
  % assemble stiffness matrix and load vector
  
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  L = assemLoad_LFE(Mesh,P7O6(),f_handle);
  
  cA = assemMat_LFE(cMesh,@STIMA_Lapl_LFE);
  
  % construct prolongation matrix
  
  P = get_PMat_LFE(cMesh,Mesh);
  
  % incorporate boundary conditions
  
  [U_bd,FDofs] = assemDir_LFE(Mesh,-1,gd_handle);
  L = L - A*U_bd;
  A = A(FDofs,FDofs);
  L = L(FDofs);
  
  [cU_bd,cFDofs] = assemDir_LFE(cMesh,-1,gd_handle);
  cA = cA(cFDofs,cFDofs);
  
  P = P(FDofs,cFDofs);
  
  % assemble matrices for scalar products
  
  H1 = A;
  L2 = assemMat_LFE(Mesh,@MASS_LFE);
  L2 = L2(FDofs,FDofs);
  
  % calculate exact solution
  
  U_ex = A\L;
  
  % define oscillatory initial guess

  u0 = @(x) cos(0.5*pi*x(:,1)).*cos(0.5*pi*x(:,2))+sin(42*pi*x(:,1)).*sin(57*pi*x(:,2))+x(:,1).^2+x(:,2).^2;
  U = u0(Mesh.Coordinates(FDofs,:));
%   U = zeros(length(FDofs),1);
  
  % calculate initial errors
  
  Err = U-U_ex;
  err_l2(1) = sqrt(Err'*L2*Err);
  err_h1(1) = sqrt(Err'*H1*Err);
  
  % do m smoothing steps
  
  for i=2:m+1

    U = smoother(U,A,L);
    
    Err = U-U_ex;
    err_l2(i) = sqrt(Err'*L2*Err);
    err_h1(i) = sqrt(Err'*H1*Err);
  end
  
  % do coarse grid correction
  
  res = P.'*(L-A*U);
  cor = cA\res;
  U = U + P*cor;
  
  Err = U-U_ex;
  err_l2(m+2) = sqrt(Err'*L2*Err);
  err_h1(m+2) = sqrt(Err'*H1*Err);
  
  % do m smoothing steps
  
  for i=m+3:2*m+2

    U = smoother(U,A,L);
    
    Err = U-U_ex;
    err_l2(i) = sqrt(Err'*L2*Err);
    err_h1(i) = sqrt(Err'*H1*Err);
  end
  
  % do coarse grid correction
  
  res = P.'*(L-A*U);
  cor = cA\res;
  U = U + P*cor;
  
  Err = U-U_ex;
  err_l2(2*m+3) = sqrt(Err'*L2*Err);
  err_h1(2*m+3) = sqrt(Err'*H1*Err);
  
  % do m smoothing steps
  
  for i=2*m+4:3*m+3

    U = smoother(U,A,L);
    
    Err = U-U_ex;
    err_l2(i) = sqrt(Err'*L2*Err);
    err_h1(i) = sqrt(Err'*H1*Err);
  end
  
  % plot errors
  
  figure;
  semilogy(0:3*m+2,err_l2/err_l2(1),'-o',0:3*m+2,err_h1/err_h1(1),'-^');
  grid on;
  set(gca,'XTick',0.5:7.5,'XTickLabel',{'GS';'GS';'CGC';'GS';'GS';'CGC';'GS';'GS'});
  legend('L^2 Error','H_0^1 Error','Location','NorthEast');
%   xlabel('\bf number of iterations');
  ylabel('\bf relative error');
  title('\bf Convergence of Gauss-Seidel Smoother with CGC');
  
return