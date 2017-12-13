% Compare smoothers for AMG

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-12;                                   % Stopping criterion
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
  
  % Define AMG options with various smoothers
  
  opt_gs = AMGDefaultOptions;
  opt_gs.pre.type = 'Gauss-Seidel CF';
  
  opt_sgs = AMGDefaultOptions;
  opt_sgs.pre.type = 'Gauss-Seidel CF';
  opt_sgs.post.type = 'Gauss-Seidel upper CF';
  opt_sgs.post.pre = 0;
  
  opt_jac = AMGDefaultOptions;
  opt_jac.pre.type = 'Jacobi';
  opt_jac.pre.its = 6;
  opt_jac.post.its = 6;
  
  opt_ainv = AMGDefaultOptions;
  opt_ainv.pre.type = 'ainv';
  opt_ainv.pre.ainv_shift = 0;
  
  % Run algebraic multigrid solver
  
  U0 = zeros(size(b));
  
  [U(FDofs),flag_gs,relres_gs,iter_gs,resvec_gs,times] = amg_solve(U0,b,A,opt_gs,TOL,MAXIT);
  time_gs = times.V_cycles.elapsed;
  
  [U(FDofs),flag_sgs,relres_sgs,iter_sgs,resvec_sgs,times] = amg_solve(U0,b,A,opt_sgs,TOL,MAXIT);
  time_sgs = times.V_cycles.elapsed;
  
  [U(FDofs),flag_jac,relres_jac,iter_jac,resvec_jac,times] = amg_solve(U0,b,A,opt_jac,TOL,MAXIT);
  time_jac = times.V_cycles.elapsed;
  
  [U(FDofs),flag_ainv,relres_ainv,iter_ainv,resvec_ainv,times] = amg_solve(U0,b,A,opt_ainv,TOL,MAXIT);
  time_ainv = times.V_cycles.elapsed;
  
  % Plot Errors
  
  figure;
  semilogy(1:iter_gs,resvec_gs,...
    1:iter_sgs,resvec_sgs,...
    1:iter_jac,resvec_jac,...
    1:iter_ainv,resvec_ainv);
  xlabel('\bf iteration');
  ylabel('\bf residual');
  title('\bf Convergence of Algebraic Multigrid');
  grid on;
  iter = max([iter_gs,iter_sgs,iter_jac,iter_ainv]);
  set(gca,'XTick',1:iter,'XLim',[1,iter],'YLim',[1e-12,1]);
  legend(sprintf('Gauss-Seidel (%3.2f s)',time_gs),...
    sprintf('GS / GS'' (%3.2f s)',time_sgs),...
    sprintf('Jacobi, 3m (%3.2f s)',time_jac),...
    sprintf('Ainv (%3.2f s)',time_ainv),...
    'Location','NorthEast');
  
  % Clear memory
  
  clear all;