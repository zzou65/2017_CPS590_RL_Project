% Test script for multigrid solver.

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  REF = [3,7];                                   % refinements of initial mesh
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-6;                                    % Stopping criterion
  MAXIT = 15;                                    % Maximum number of iterations
  CYC = 1;                                       % V or W multigrid cycles
  M = [1,1];                                     % Smoothing steps 
  SMOOTHER = @gs_smooth;                         % Multrigrid smoother
  
  % Generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',REF);
  mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE);
  mg_data = mg_smooth(mg_data,'m',M,'cyc',CYC,'smoother',SMOOTHER);
  mg_data = mg_error(mg_data,'energy',false,'l2',true,'rel',true);
  
  hmg_data = mg_smooth(mg_data,'type','hier');
  
  % Generate algebraic multigrid options
  
  AMGopt = AMGDefaultOptions;
  AMGopt.pre.its = M(1);
  AMGopt.post.its = M(end);
  
  A = mg_data{end}.A;
  b = mg_data{end}.b;
  
  % Run multigrid solvers
  
  [x_mg,conv_mg] = mg(mg_data,TOL,MAXIT);
  resvec_mg = conv_mg.error.l2_iter/conv_mg.error.l2_iter(1);
  [x_hmg,conv_hmg] = mg(hmg_data,TOL,MAXIT);
  resvec_hmg = conv_hmg.error.l2_iter/conv_hmg.error.l2_iter(1);
  [x_amg,flag_amg,relres_amg,iter_amg,resvec_amg] = amg_solve(zeros(size(b)),b,A,AMGopt,TOL,MAXIT);
  resvec_amg = resvec_amg/resvec_amg(1);
  
  % Plot errors
  
  figure;
  semilogy(1:conv_mg.iter,resvec_mg,'-o',...
    1:conv_hmg.iter,resvec_hmg,'-x',...
    1:iter_amg,resvec_amg,'-+');
  xlabel('\bf iteration');
  ylabel('\bf residual');
  title('\bf Convergence of Multigrid Methds');
  legend('multigrid','hierarchical multigrid','algebraic multigrid','Location','NorthEast');
  grid on;
  
  % Clear memory
  
  clear all;
  