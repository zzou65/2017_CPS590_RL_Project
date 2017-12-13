function [] = complexity()
% complexity of nested multigrid
%
%   plots the time required by nested multigrid to solve a sample problem
%   vs. the number of degrees of freedom on the fines mesh for various
%   amounts of refinement and the least-squares line through these points.
%   The slope of this line is approximately one, so the algorithm is
%   optimal.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize constants
 
  ref = [3 8];                                   % refinements of mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  du_ex = @(x,varargin) 2*x;                     % derivative of exact solution
  smoother = @gs_smooth;                         % multigrid smoother
  tol_c = 1e-2;                                  % Stopping criterion
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
  
  % determine number of degrees of freedom on all levels
  
  dofs = zeros(1,LVL-1);
  for lvl=2:LVL
    dofs(lvl-1) = mg_data{lvl}.n.free;
  end
  
  % time nested multigrid for various levels of refinement
  
  time = zeros(1,LVL-1);
  for lvl=2:LVL
    
    % calculate exact discrete solution and discretization error
    u = mg_data{lvl}.u_bd;
    u(mg_data{lvl}.dofs) = mg_data{lvl}.A\mg_data{lvl}.b;
    disc_err = H1SErr_LFE(mg_data{lvl}.mesh,u,P7O6(),du_ex);
    tol = tol_c*disc_err;
    
    % determine reasonable number of times to calculate solution
    %   (timer is imprecise for short times)
    num = ceil(1e4*exp(-log(dofs(lvl-1))));
    
    % start timer
    t = cputime;
    
    % run nested multigrid solver
    for i=1:num
      nmg_solve({mg_data{1:lvl}},tol,maxit);
    end
    
    % stop timer
    time(lvl-1) = (cputime - t)/num;
    
  end
  
  % plot time against number of degrees of freedom
  
  figure;
  plot(dofs,time,'bo');
  hold on;
  set(gca,'XScale','log','YScale','log');
  grid('on');

  title('{\bf Complexity of Nested Multigrid}');
  xlabel('{\bf degrees of freedom}');
  ylabel('{\bf time}')

  p = polyfit(log(dofs),log(time),1);
  x = get(gca,'XLim');
  y = get(gca,'YLim');
  plot(x,exp(polyval(p,log(x))));
  set(gca,'YLim',y);
  add_Slope(gca,'NorthWest',p(1));
  
return