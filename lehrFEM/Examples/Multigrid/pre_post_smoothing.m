function [] = pre_post_smoothing()
% comparison of presmoothing and postsmoothing
%
%   This code solves a sample problem with the multigrid methid using a
%   total of five smoothing steps, divided in various ways between
%   presmoothing and postsmoothing.  The error in the energy norm is
%   plotted after every V-cycle.
%
%   The convergence of the multigrid iteration is fastest when the
%   smoothing iterations are diveded evenly between presmoothing and
%   postsmoothing, possibly with slightly more presmoothing than
%   postsmoothing.  This tendency, however, is switched if the order in the
%   Gauss-Seidel smoothing iteration is reversed.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  ref = [0 7];                                   % refinements of mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  smoother = @gs_smooth;                         % multigrid smoother
  tol = 1e-6;                                    % Stopping criterion
  maxit = 100;                                   % Maximum number of iterations
  cyc = 1;                                       % V or W multigrid cycles
  m = 5;                                         % total number of smoothing steps
  
  % generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',ref);
  mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
  mg_data = mg_smooth(mg_data,'cyc',cyc,'smoother',smoother);
%   mg_data = mg_smooth(mg_data,'cyc',cyc,'smoother',smoother,'per','FC'); % use this line to switch order of smoothing
  mg_data = mg_error(mg_data,'iter',false,'exact',true,'rel',false);
  LVL = length(mg_data);
  
  % initialize data
  
  numit = nan(1,m+1);
  err = cell(1,m+1);
  
  x0 = zeros(size(mg_data{LVL}.b));
  err0 = mg_data{LVL}.error.energy_exact(x0);
  
  % run multigrid solver for various divisions of m into presmoothing and
  % postsmoothing
  
  for m1=0:m
    mg_data = mg_smooth(mg_data,'m',[m1,m-m1]);
    [x,conv] = mg(x0,mg_data,tol,maxit);
    numit(m1+1) = conv.iter;
    err{m1+1} = [err0,conv.error.energy_exact];
  end
  
  % plot error versus number of v-cycles
  
  figure;
  lgd = cell(1,m+1);
  for i=1:m+1
    semilogy(0:numit(i),err{i},'-o');
    hold all;
    lgd{i} = sprintf('m_{pre} = %.0f, m_{post} = %.0f',i-1,m+1-i);
  end
  legend(lgd{:},'Location','NorthEast');
  grid on;
  set(gca,'YLim',[tol,err0]);
  xlabel('\bf V-cycles');
  ylabel('\bf error in energy norm');
  title('\bf Multigrid with Gauss-Seidel Presmoothing and Postsmoothing');

return