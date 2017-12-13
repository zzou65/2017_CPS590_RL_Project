function [] = mg_full_hier()
% compare (unnested) full multigrid and hierarchical multigrid
%
%   This code solves a sample problem with full and hierarchical multigrid
%   solvers and displays the convergence data in two figures.  In the
%   first, the energy-norm error is plotted after every V-cycle against the
%   number of 1-dimensional smoothing steps done up to that iteration. The
%   second figure shows the total amount of time spent on smoothing and on
%   intergrid transfer by the two methods.
%
%   Both figures clearly show that hierarchical multigrid is much less
%   efficient than full multigrid.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize constants
 
  ref = [0 7];                                   % refinements of mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  smoother = @gs_smooth;                         % multigrid smoother
  tol = 1e-3;                                    % Stopping criterion
  maxit = 100;                                   % Maximum number of iterations
  cyc = 1;                                       % V or W multigrid cycles
  m = [1 1];                                     % number of smoothing steps

  % generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',ref);
  mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
  mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother);
  mg_data = mg_error(mg_data,'iter',false,'exact',true,'rel',false);
  LVL = length(mg_data);
  
  hmg_data = mg_smooth(mg_data,'type','hier');
  
  % determine cost of a v-cycle
  
  numssc0_mg = 0;
  numssc0_hmg = 0;
  for i=2:LVL
    numssc0_mg = numssc0_mg + mg_data{i}.n.smooth;
    numssc0_hmg = numssc0_hmg + hmg_data{i}.n.smooth;
  end
  
  % initialize data
  
  x0 = zeros(size(mg_data{LVL}.b));
  err0 = mg_data{LVL}.error.energy_exact(x0);
  
  timetransf_mg = 0;
  timessc_mg = 0;
  timetransf_hmg = 0;
  timessc_hmg = 0;
  
  runs = 10;
  
  % run full multigrid solver
  
  for j=1:runs
    [x,conv_mg] = mg(x0,mg_data,tol,maxit);
    timetransf_mg = timetransf_mg + conv_mg.time.transfer/runs;
    timessc_mg = timessc_mg + conv_mg.time.smooth/runs;
  end
  numssc_mg = numssc0_mg*(0:conv_mg.iter);
  err_mg = [err0,conv_mg.error.energy_exact];
  
  % run hierarchical multigrid solver
  
  for j=1:runs
    [x,conv_hmg] = mg(x0,hmg_data,tol,maxit);
    timetransf_hmg = timetransf_hmg + conv_hmg.time.transfer/runs;
    timessc_hmg = timessc_hmg + conv_hmg.time.smooth/runs;
  end
  numssc_hmg = numssc0_hmg*(0:conv_hmg.iter);
  err_hmg = [err0,conv_hmg.error.energy_exact];
  
  % plot error versus number of 1-d subspace corrections
  
  figure;
  semilogy(numssc_mg,err_mg,'-o',numssc_hmg,err_hmg,'-v');
  legend('full multigrid','hierarchical multigrid','Location','NorthEast');
  grid on;
  set(gca,'YLim',[tol,err0]);
  xlabel('\bf 1-dimensional subspace corrections');
  ylabel('\bf error in energy norm');
  title('\bf Full and Hierarchical Multigrid');
  
  % plot time for subspace corrections and intergrid transfer
  
  figure;
  bar([1 2],[timessc_mg,timetransf_mg;timessc_hmg,timetransf_hmg]);
  legend('Subspace Correction','Intergrid Transfer','Location','NorthWest');
  set(gca,'XTick',[1 2],'XTickLabel',{'Full Multigrid';'Hierarchical Multigrid'},'YGrid','on');
  ylabel('\bf time to reach tolerance 10^{-3} in energy norm');
  title('\bf Full and Hierarchical Multigrid');
  
return