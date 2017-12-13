% Test script for algebraic multigrid solver

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-6;                                    % Stopping criterion
  MAXIT = 50;                                    % Maximum number of iterations
 
  % Generate coarse mesh

  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
  % Refine mesh
  
  for i = 1:5
    Mesh = refine_REG(Mesh);
  end
  
  % Construct stiffness matrix and load vector
  
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  b = assemLoad_LFE(Mesh,P3O3(),F_HANDLE);
  
  % Incoporate Dirichlet boundary conditions
  
  [U,FDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
  b = b - A*U;
  A = A(FDofs,FDofs);
  b = b(FDofs);
  
  % Define AMG options
  
  AMGOptions = AMGDefaultOptions;
  
  AMGOptions.pre.its = 1;
  AMGOptions.post.its = 1;
  AMGOptions.CF.theta = 0.2;
  AMGOptions.P.method64 = 'original';
  
  % Run algebraic multigrid solver
  
  t = cputime;
  U0 = zeros(size(b));
  [U(FDofs),flag,relres,iter,resvec,times] = amg_solve(U0,b,A,AMGOptions,TOL,MAXIT);
  time = cputime - t;
  fprintf('Total Time : %g\n',time);
  fprintf('Time for V-cycles : %g\n',times.V_cycles.elapsed);
  
  % Plot errors
  
  figure;
  semilogy(1:iter,resvec);
  xlabel('\bf iteration');
  ylabel('\bf residual');
  title('\bf Convergence of Algebraic Multigrid');
  grid on;
  set(gca,'XTick',1:iter,'XLim',[1,iter]);
  
  p = polyfit(1:iter,log(resvec),1);
  add_Slope(gca,'SouthWest',p(1));
  
  % Clear memory
  
  clear all;
  