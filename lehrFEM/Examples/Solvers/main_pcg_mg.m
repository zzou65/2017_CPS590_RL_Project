% Run script for CG solver with multigrid preconditioner.

%   Copyright 2005-2007 Patrick Meury & Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  REF = [2 7];                                   % Number of red refinements
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-12;                                   % Stopping criterion
  MAXIT = 100;                                   % Maximum number of iterations
  IT = 1;                                        % Multigrid steps
  CYC = 1;                                       % V or W multigrid cycles
  M = 1;                                         % Initial relaxation steps 

  % Construct multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',REF);
  mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE);
  mg_data = mg_smooth(mg_data,'m',M,'cyc',CYC,...
    'smoother',@gs_smooth,'postsmoother',@gst_smooth);
  
  A = mg_data{end}.A;
  L = mg_data{end}.b;
  U = mg_data{end}.u_bd;
  
  % Run preconditioned CG solver
 
  t = cputime;
  U0 = zeros(size(L));
  [U(mg_data{end}.dofs),flag,relres,iter,resvec] = ...
    pcg_solve(U0,A,L,TOL,MAXIT,@mg,IT,mg_data);
  fprintf('Runtime of preconditioned CG solver [s]  :  %f\n',cputime-t);
  
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
  