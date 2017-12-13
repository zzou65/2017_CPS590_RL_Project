% Run script for multigrid solver with nested iterations.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % Initialize constants
  
  REF = [2 8];                                   % Number of red refinements
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-12;                                   % Stopping criterion
  MAXIT = 100;                                   % Maximum number of iterations
  CYC = 1;                                       % V or W multigrid cycles
  M = 1;                                         % Relaxation steps
  
  % Generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',REF);
  mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE);
  mg_data = mg_smooth(mg_data,'m',M,'cyc',CYC,'smoother',@gs_smooth);
  mg_data = mg_error(mg_data,'energy',false,'eucl',true);
  LVL = length(mg_data);
  
  % calculate tolerances for multigird
  % (so that the stopping criteria are the same as in other methods)
    
  tol = TOL*4.^(LVL-2:-1:0);
  for i=1:LVL-1
    nrmb = norm(mg_data{i+1}.b);
    if(nrmb == 0)
      nrmb = 1;
    end
    tol(i) = tol(i)*nrmb;
  end
  
  % Solve with nested multigrid
  
  U = mg_data{end}.u_bd;
	[U(mg_data{end}.dofs),conv] = nmg_solve(mg_data,tol,MAXIT);
  
  % display convergence information
  
  total = conv.time.coarse_solve0;
  for i = 1:LVL-1
  
    t = conv.time.smooth(i)+conv.time.transfer(i)+conv.time.coarse_solve(i);
    total = total + t;
    
    if(~conv.flag(i))
      error('Nested multigrid solver failed to converge.');  
    else
      fprintf('\n');
      fprintf('Level %2d\n',i+1);
      fprintf('  Runtime [s]           :  %.1e\n',t);
      fprintf('  Tolerance             :  %.1e\n',conv.tol(i));
      fprintf('  Multigrid iterations  :  %d\n',conv.iter(i));
      fprintf('\n');
    end
    
  end
  
  fprintf('Total runtime [s]  :  %.1e\n',total);
  
  % Clear memory
  
  clear all;
  