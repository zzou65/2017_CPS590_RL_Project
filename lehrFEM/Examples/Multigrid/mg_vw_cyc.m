function [] = mg_vw_cyc()
% compare multigrid v-cycles and w-cycles
%
%   This code solves a sample problem with multigrid solvers using V-cycles
%   and W-cycles and displays the convergence data in two figures.  In the
%   first, the energy-norm error is plotted after every iteration against
%   the number of 1-dimensional smoothing steps done up to that iteration.
%   The second figure shows the total amount of time spent on smoothing and
%   on intergrid transfer by the two methods.
%
%   In this experiment, W-cycles are somewhat more expensive than V-cycles.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize constants
 
  ref = [0 7];                                   % refinements of mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  smoother = @gs_smooth;                         % multigrid smoother
  tol = 1e-6;                                    % Stopping criterion
  maxit = 100;                                   % Maximum number of iterations
  m = [1 1];                                     % number of smoothing steps

  % generate multigrid data structure for v-cycles
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data_v = mg_mesh('mesh',CMesh,'ref',ref);
  mg_data_v = mg_stima(mg_data_v,'f',f_handle,'gd',gd_handle);
  mg_data_v = mg_smooth(mg_data_v,'m',m,'cyc',1,'smoother',smoother);
  mg_data_v = mg_error(mg_data_v,'iter',false,'exact',true,'rel',false);
  LVL = length(mg_data_v);
  
  % adapt data structure for w-cycles
  
  mg_data_w = mg_smooth(mg_data_v,'cyc',2);
  
  % determine cost of a v-cycle
  
  numssc0_v = 0;
  for i=2:LVL
    numssc0_v = numssc0_v + mg_data_v{i}.n.smooth;
  end
  
  % determine cost of a w-cycle
  
  numssc0_w = 0;
  for i=2:LVL
    numssc0_w = 2*numssc0_w + mg_data_w{i}.n.smooth;
  end
  
  % initialize data
  
  x0 = zeros(size(mg_data_v{LVL}.b));
  err0 = mg_data_v{LVL}.error.energy_exact(x0);
  
  runs = 10;
  
  timetransf_v = 0;
  timessc_v = 0;
  timetransf_w = 0;
  timessc_w = 0;
  
  % run multigrid solver with v-cycles
  
  for j=1:runs
    [x,conv_v] = mg(x0,mg_data_v,tol,maxit);
    timetransf_v = timetransf_v + conv_v.time.transfer/runs;
    timessc_v = timessc_v + conv_v.time.smooth/runs;
  end
  numssc_v = numssc0_v*(0:conv_v.iter);
  err_v = [err0,conv_v.error.energy_exact];
  
  % run multigrid solver with v-vycles
  
  for j=1:runs
    [x,conv_w] = mg(x0,mg_data_w,tol,maxit);
    timetransf_w = timetransf_w + conv_w.time.transfer/runs;
    timessc_w = timessc_w + conv_w.time.smooth/runs;
  end
  numssc_w = numssc0_w*(0:conv_w.iter);
  err_w = [err0,conv_w.error.energy_exact];
  
  % plot error versus number of multigrid iterations (v-cycles or w-cycles)
  
  figure;
  semilogy(0:conv_v.iter,err_v,'-v',0:conv_w.iter,err_w,'-s');
  legend('V-cycles','W-cycles','Location','NorthEast');
  grid on;
  set(gca,'YLim',[tol,err0]);
  xlabel('\bf multigrid iterations');
  ylabel('\bf error in energy norm');
  title('\bf Multigrid V-Cycles and W-Cycles');

  % plot error versus number of 1-d subspace corrections
  
  figure;
  semilogy(numssc_v,err_v,'-v',numssc_w,err_w,'-s');
  legend('V-cycles','W-cycles','Location','NorthEast');
  grid on;
  set(gca,'YLim',[tol,err0]);
  xlabel('\bf 1-dimensional subspace corrections');
  ylabel('\bf error in energy norm');
  title('\bf Multigrid V-Cycles and W-Cycles');
  
  % plot time for subspace corrections and intergrid transfer
  
  figure;
  bar([1 2],[timessc_v,timetransf_v;timessc_w,timetransf_w]);
  legend('Subspace Correction','Intergrid Transfer','Location','NorthWest');
  set(gca,'XTick',[1 2],'XTickLabel',{'V-cycles';'W-cycles'},'YGrid','on');
  ylabel('\bf time to reach tolerance 10^{-6} in energy norm');
  title('\bf Multigrid V-Cycles and W-Cycles');
  
return  