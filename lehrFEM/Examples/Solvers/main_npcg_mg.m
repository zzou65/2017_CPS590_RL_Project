% Run script for nested CG solver with multigrid preconditioner.

%   Copyright 2005-2007 Patrick Meury & Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  IREFS = 2;                                     % Inital red refinements (2 or greater)
  NREFS = 6;                                     % Number of red refinements
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-12;                                   % Stopping criterion
  MAXIT = 100;                                   % Maximum number of iterations
  IT = 1;                                        % Multigrid steps
  CYC = 1;                                       % V or W multigrid cycles
  M = 1;                                         % Initial relaxation steps
  
  % Set up prescribed tolerance
  
  tol = TOL*4^(IREFS+NREFS-1);
  
  % Generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',[IREFS,IREFS+NREFS]);
  mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE);
  mg_data = mg_smooth(mg_data,'m',M,'cyc',CYC,...
    'smoother',@gs_smooth,'postsmoother',@gst_smooth);
  
  % Solve on coarsest grid
  
  total = 0;
  t0 = cputime();
  
  U0 = mg_data{1}.A\mg_data{1}.b;
  
  dt = cputime() - t0;
  total = total + dt;
  
  for lvl = 2:NREFS+1
    
    % Prolong coarse grid solution to fine grid
    
    U = mg_data{lvl-1}.u_bd;
    U(mg_data{lvl-1}.dofs) = U0;
    U0 = mg_data{lvl}.P_full(mg_data{lvl}.dofs,:)*U;
    
    % Compute solution on current level
    
    A = mg_data{lvl}.A;
    L = mg_data{lvl}.b;
    
    t0 = cputime();
    [U0,flag,relres,iter,resvec] = ...
      pcg_solve(U0,A,L,TOL,MAXIT,@mg,IT,{mg_data{1:lvl}});
    dt = cputime() - t0;
    total = total + dt;
    if(flag == 0)
      error('PCG solver failed to converge.');  
    else
      fprintf('\n');
      fprintf('Level %2d\n',lvl);
      fprintf('  Runtime [s]         :  %.1e\n',dt);
      fprintf('  Relative tolerance  :  %.1e\n',tol);
      fprintf('  PCG iterations      :  %d\n',iter);
      fprintf('\n');
    end
    
  end
    
  fprintf('Total runtime [s]  :  %.1e\n',total);
 
  % Generate figure
  
  fig = figure('Name','Preconditioned CG solver');
  plot(resvec,'rx');
  if(CYC == 1)
    title(['{\bf Preconditioned CG solver (' int2str(IT) '-step V(' ...
           int2str(M) ',' int2str(M) ') preconditioner)}']);
  else
    title(['{\bf Preconditioned CG solver (' int2str(IT) '-step W(' ...
           int2str(M) ',' int2str(M) ') preconditioner)}']);  
  end
  xlabel('{\bf Iteration number}');
  ylabel('{\bf Relative residual [log]}');
  set(gca,'YScale','log');
  
  % Clear memory
  
  clear all;
  