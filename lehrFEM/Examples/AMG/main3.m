% discontinuous coefficients
%
%   Compares multigrid and algebraic multigrid for discontinuous
%   conductivity

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  REF = [2,6];                                   % refinements of initial mesh
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-6;                                    % Stopping criterion
  MAXIT = 15;                                    % Maximum number of iterations
  CYC = 1;                                       % V or W multigrid cycles
  M = [1,1];                                     % Smoothing steps 
  NRM = 0;
  SMOOTHER = @gs_smooth;                         % Multrigrid smoother
  
  % Define conductivity
  
  epsilon = 1e-5;
  C_HANDLE = @(x,varargin) 1*(x(:,1)<0.5)...
                        + 0.5*(1+epsilon)*(x(:,1)==0.5)...
                        + epsilon*(x(:,1)>0.5);
  
  % Generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  
  mg_data = mg_mesh('mesh',CMesh,'ref',REF);
  mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE,'c',C_HANDLE);
  mg_data = mg_smooth(mg_data,'m',M,'cyc',CYC,'smoother',SMOOTHER);
  mg_data = mg_error(mg_data,'energy',false,'eucl',true);
  
  A = mg_data{end}.A;
  b = mg_data{end}.b;
  
  % Generate algebraic multigrid options
  
  AMGopt = AMGDefaultOptions;
  AMGopt.pre.its = M(1);
  AMGopt.post.its = M(end);
  AMGopt.mincoarse = mg_data{2}.n.free-1;

  % Run multigrid solvers

  nrmb = norm(mg_data{end}.b);
  if(nrmb == 0)
    nrmb = 1;
  end
  [x_mg,conv_mg] = mg(mg_data,TOL*nrmb,MAXIT);
  [x_amg,flag_amg,relres_amg,iter_amg,resvec_amg] = amg_solve(zeros(size(b)),b,A,AMGopt,TOL,MAXIT);
  
  % Plot errors
  
  figure;
  semilogy(1:conv_mg.iter,conv_mg.error.eucl_iter/nrmb,'-o',... %1:iter_mg,resvec_mg,'-o',...
    1:iter_amg,resvec_amg,'-+');
  xlabel('\bf iteration');
  ylabel('\bf residual');
  title(sprintf('%s Convergence of Multigrid Methds for Discontinuous Coefficients (%s = %g)','\bf','\epsilon',epsilon));
  legend('multigrid','algebraic multigrid','Location','NorthEast');
  grid on;
  
  % Clear memory
  
  clear all;
  