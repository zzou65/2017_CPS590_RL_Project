function [] = presmoothing_its()
% effect of number of presmoothing steps
%
%   This code solves a sample problem using the multigrid method with
%   various amounts of Gauss-Seidel presmoothing and displays the
%   convergence data in two figures.  In the first, the energy-norm error
%   is plotted after every V-cycle against the number of 1-dimensional
%   smoothing steps done up to that iteration. The second figure shows the
%   total amount of time spent on smoothing and on intergrid transfer by
%   the multigrid methods.
%
%   The results suggest that it is most efficient to use only one smoothing
%   step per iteration for this problem.

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
  
  % define values for m(1)
  
  m = 1:5;
  num = length(m);
  
  % generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
%   CMesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',ref);
  mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
  mg_data = mg_smooth(mg_data,'m',[1 0],'cyc',cyc,'smoother',smoother);
  mg_data = mg_error(mg_data,'iter',false,'exact',true,'rel',false);
  LVL = length(mg_data);
  
  % determine cost of a v-cycle
  
  numssc0 = 0;
  for i=2:LVL
    numssc0 = numssc0 + mg_data{i}.n.smooth;
  end
  
  % initialize data
  
  numssc = cell(1,num);
  err = cell(1,num);
  timessc = zeros(1,num);
  timetransf = zeros(1,num);
  
  x0 = zeros(size(mg_data{LVL}.b));
  err0 = mg_data{LVL}.error.energy_exact(x0);
  
  runs = 25;
  
  % run multigrid solver for various m
  
  for i=1:num
    mg_data = mg_smooth(mg_data,'m',[m(i) 0]);
    for j=1:runs
      [x,conv] = mg(x0,mg_data,tol,maxit);
      
      timetransf(i) = timetransf(i) + conv.time.transfer/runs;
      timessc(i) = timessc(i) + conv.time.smooth/runs;
    end
    numssc{i} = numssc0*m(i)*(0:conv.iter);
    err{i} = [err0,conv.error.energy_exact];
  end
  
  % plot error versus number of 1-d subspace corrections
  
  figure;
  lgd = cell(1,num);
  for i=1:num
    semilogy(numssc{i},err{i},'-o');
    hold all;
    lgd{i} = sprintf('m = %.0f',m(i));
  end
  legend(lgd{:},'Location','NorthEast');
  grid on;
  set(gca,'YLim',[tol,err0],'XLim',[0,7e5]);
  xlabel('\bf 1-dimensional subspace corrections');
  ylabel('\bf error in energy norm');
  title('\bf Multigrid with Gauss-Seidel Presmoothing');
  
  % plot time for subspace corrections and intergrid transfer vs m
  
  figure;
  bar(m,[timessc;timetransf]');
  legend('Subspace Correction','Intergrid Transfer','Location','NorthWest');
  set(gca,'XTick',m,'YGrid','on');
  xlabel('\bf presmoothing steps m');
  ylabel('\bf time to reach tolerance 10^{-6} in energy norm');
  title('\bf Multigrid with Gauss-Seidel Presmoothing');
  
return