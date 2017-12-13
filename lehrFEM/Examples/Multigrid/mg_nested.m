function [] = mg_nested()
% compare nested and unnested multigrid iterations
%
%   This code solves a sample problem with nested and unnested multigrid
%   solvers, using a random initial guess for the unnested method, and
%   displays the convergence data in two figures.  In the first, the
%   energy-norm error is plotted after every V-cycle against the number of
%   1-dimensional smoothing steps done up to that iteration. The second
%   figure shows the total amount of time spent on smoothing and on
%   intergrid transfer by the two methods.  For nested multigrid, the error
%   is only plotted once the iteration reaches the finest grid.
%
%   This experiment shows that nested multigrid is superior to unnested
%   multigrid with a random initial guess, using just two thirds of the
%   time of the latter.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % initialize constants
 
  ref = [0 7];                                   % refinements of mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
%   f_handle = @f_LShap;                           % Right hand-side source term
%   gd_handle = @g_D_LShap;                        % Dirichlet boundary data
  smoother = @gs_smooth;                         % multigrid smoother
  tol = 1e-6;                                    % Stopping criterion
  maxit = 100;                                   % Maximum number of iterations
  cyc = 1;                                       % V or W multigrid cycles
  m = [1 1];                                     % number of smoothing steps

  % generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
%   CMesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',ref);
  mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
  mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother);
  mg_data = mg_error(mg_data,'iter',false,'exact',true,'rel',false);
  LVL = length(mg_data);
  
  % determine cost of a v-cycle on each level
  
  numssc0 = zeros(1,LVL-1);
  for i=2:LVL
    numssc0(i-1:LVL-1) = numssc0(i-1:LVL-1) + size(mg_data{i}.A,1);
  end
  
  % initialize data
  
  rand('state',3);
  x0 = rand(size(mg_data{LVL}.b));
  err0 = mg_data{LVL}.error.energy_exact(x0);
  
  runs = 10;
  
  timetransf = 0;
  timessc = 0;
  timetransf_n = 0;
  timessc_n = 0;
  
  % run multigrid solver
  
  for j=1:runs
    [x,conv] = mg(x0,mg_data,tol,maxit);
    timetransf = timetransf + conv.time.transfer/runs;
    timessc = timessc + conv.time.smooth/runs;
  end
  numssc = numssc0(LVL-1)*(m(1)+m(end))*(0:conv.iter);
  err = [err0,conv.error.energy_exact];
  
  % run nested multigrid solver
  
  for j=1:runs
    [x,conv] = nmg_solve(mg_data,tol,maxit);
    timetransf_n = timetransf_n + sum(conv.time.transfer)/runs;
    timessc_n = timessc_n + sum(conv.time.smooth)/runs;
  end
  numssc_c = (m(1)+m(end))*numssc0(1:LVL-2)*conv.iter(1:LVL-2)';
  numssc_n = numssc_c+(m(1)+m(end))*numssc0(LVL-1)*(1:conv.iter(LVL-1));
  err_n = conv.error{LVL-1}.energy_exact;


  % plot error versus number of 1-d subspace corrections
  
  figure;
  semilogy(numssc,err,'-o',numssc_n,err_n,'-^');
  legend('Multigrid','Nested Multigrid','Location','NorthEast');
  grid on;
  set(gca,'YLim',[tol,err0]);
  xlabel('\bf 1-dimensional subspace corrections');
  ylabel('\bf error in energy norm');
  title('\bf Nested and Unnested Multigrid');
  
  % plot time for subspace corrections and intergrid transfer
  
  figure;
  bar([1 2],[timessc,timetransf;timessc_n,timetransf_n]);
  legend('Subspace Correction','Intergrid Transfer','Location','NorthEast');
  set(gca,'XTick',[1 2],'XTickLabel',{'Multigrid';'Nested Multigrid'},'YGrid','on');
  ylabel('\bf time to reach tolerance 10^{-6} in energy norm');
  title('\bf Nested and Unnested Multigrid');
  
return